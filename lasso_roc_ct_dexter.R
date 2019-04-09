library(glmnet)
library(stringr)
rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)

ct = commandArgs(trailingOnly=TRUE)[1]
nb_nucl = commandArgs(trailingOnly=TRUE)[2]
load("Data/unique_annot_ct.RData")

ct = "CNhs11252"
#nb_nucl = 2


# Variables de composition nucléotidiques pour les enhancers
dexter_vars_train <- read.table(paste0("Dexter_results/", ct,"/models/all_domains.dat_training_set.log.matrix"), h = T, sep = ' ')
dexter_vars_test <- read.table(paste0("Dexter_results/", ct,"/models/all_domains.dat_testing_set.log.matrix"), h = T, sep = ' ')

dexter_vars <- rbind.data.frame(dexter_vars_train, dexter_vars_test, make.row.names = TRUE )

row.names(dexter_vars) = dexter_vars$sequence
dexter_vars = dexter_vars[ ,(names(dexter_vars) != 'sequence' & names(dexter_vars)!= "expression")]


#Jeu de données comprenant les paires enhancers-TSS positives et négatives
pairs = read.table(paste0("EG_pairs/EG_pairs_", ct, ".txt"), sep = ';', h = T)
cols = c("pairs", paste0("E_", names(dexter_vars)),"chr", "position", "distance", "interaction")
data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
colnames(data) <- cols
data$interaction = pairs$Interaction
data$distance = pairs$Distance
data$position = pairs$Enhancer_center
data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]

#Pour chaque paire, leur attribuer les variables de taux de nucléotides
for(i in 1:length(pairs$Enhancer)){
  data$pairs[i] = paste0(pairs$Enhancer[i],"/", pairs$Gene[i])
  data[i, paste0("E_", names(dexter_vars))] = dexter_vars[pairs$Enhancer[i],]
}
#order by chromosomal position
data <- data[order(data$pairs),]
rownames(data) <- data$pairs
data <- na.omit(data)

#Separate training set and evaluation set
data_train <- data[1:round(dim(data)[1]*4/5),]
data_test <- data[which(!rownames(data) %in% rownames(data_train)),]


#defining variables
variables_ep = names(data)[which(names(data)!="interaction" & names(data)!="pairs" & names(data)!="chr"
                                 & names(data)!="position")]

#variables_e = variables_ep[which(!variables_ep %in% paste0("P_", row.names(quadri_proms)))]

y <- as.matrix(data_train$interaction)

#model with enhancers and promoter variables
x_ep <- as.matrix(data_train[, variables_ep])
glm_ep <- cv.glmnet(x_ep, y, alpha = 1, family="binomial", keep = T, nfolds = 10)

#model with enhancer variables
#x_e  <- as.matrix(data_train[, variables_e])
#glm_e <- cv.glmnet(x_e, y, alpha = 1, family="binomial", keep = T, nfolds = 10)

#evaluate the models and compare them
x_test_ep <- as.matrix(data_test[, variables_ep])
#x_test_e <- as.matrix(data_test[, variables_e])

y_test <- as.matrix(data_test$interaction)

pred_ep <- as.numeric(predict(glm_ep, x_test_ep, type = "response"))
#pred_e <- as.numeric(predict(glm_e, x_test_e, type = "response"))


iLambdaMinq_ep = which(glm_ep$lambda==glm_ep$lambda.min)
#iLambdaMinq_e = which(glm_e$lambda==glm_e$lambda.min)

positivesrates = function(seuil, d, pred){
  pbin = ifelse(pred>seuil, 1, 0)
  err = d$interaction - pbin
  err2 = d$interaction + pbin
  
  TP = sum(err2==2)
  FP = sum(err == -1)
  c(TP/sum(d$interaction==1), FP/sum(d$interaction==0))
}


library(pROC) 
auc_ep <- roc(data_test$interaction, pred_ep)$auc
#auc_e <- roc(data_test$interaction, pred_e)$auc

curve_ep = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = pred_ep)
#curve_e = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = pred_e)

png(paste0('Figures/roc',ct, '_dexter.png'))
plot(curve_ep[2,], curve_ep[1,], type = "l", col = "blue", 
     main = paste("ROC", ct, " dexter AUC : ", as.character(round(auc_ep,4))), xlim = c(0,1), ylim = c(0,1),
     xlab = "FP rate", ylab = "TP rate")
#lines(curve_e[2,], curve_e[1,], type = "l", col = "red")
#legend(0.3, 0.2, legend=c(paste0("Enhancers and promoters features, auc = ", as.character(round(auc_ep,4))), 
                          #paste0("Enhancers features, auc = ", as.character(round(auc_e,4)))),
       col=c("blue", "red"), lty=1:2, cex=0.8)
abline(c(0,0),c(1,1))
dev.off()

varselq =names(which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0))
glm_ep$glmnet.fit$beta[,iLambdaMinq_ep][which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0)]

data_description <- function(data_train){
  par(mfrow = c(2,2))
  enh_train_pos <- str_split_fixed(subset(data_train, data_train$interaction == 1)$pairs, '/', 2)[,1]
  prom_train_pos <- str_split_fixed(subset(data_train, data_train$interaction == 1)$pairs, '/', 2)[,2]
  
  hist(table(enh_train_pos), breaks = 20, xlim = c(0,15),
       main = paste("# genes per enhancers in", ct), col = "darksalmon")
  hist(table(prom_train_pos), breaks = 20, xlim = c(0,15), col = "cadetblue3",
       main = "# enhancers per gene, Positives")
  
  enh_train_neg <- str_split_fixed(subset(data_train, data_train$interaction == 0)$pairs, '/', 2)[,1]
  prom_train_neg <- str_split_fixed(subset(data_train, data_train$interaction == 0)$pairs, '/', 2)[,2]
  
  hist(table(enh_train_neg), breaks = 20, xlim = c(0,15),
       main = paste("# genes per enhancers in", ct), col = "darksalmon")
  hist(table(prom_train_neg), breaks = 20, xlim = c(0,15), col = "cadetblue3",
       main = "# enhancers per gene, Negatives")
}

train_test_intersect <- function(){
  enh_train <- str_split_fixed(data_train$pairs, '/', 2)[,1]
  prom_train <- str_split_fixed(data_train$pairs, '/', 2)[,2]
  
  enh_test <- str_split_fixed(data_test$pairs, '/', 2)[,1]
  prom_test <- str_split_fixed(data_test$pairs, '/', 2)[,2]
  
  length(which(enh_test %in% enh_train))
  length(which(prom_test %in% prom_train))
  prom_test[which(prom_test %in% prom_train)]
}

#data_description(data)
