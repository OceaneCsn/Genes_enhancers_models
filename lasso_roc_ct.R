library(glmnet)
library(stringr)
rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)

ct = commandArgs(trailingOnly=TRUE)[1]
nb_nucl = commandArgs(trailingOnly=TRUE)[2]
load("Data/unique_annot_ct.RData")

ct = "CNhs11344"
nb_nucl = 2

# Variables de composition nucléotidiques pour les enhancers
quadri_enhancers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl),".txt"), sep = ';', h = T, na.strings="")
headers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl), ".txt"), sep = ';', h = F, na.strings="")[1,]
row.names(quadri_enhancers) = quadri_enhancers$X
quadri_enhancers = quadri_enhancers[ ,(names(quadri_enhancers) != 'X')]
colnames(quadri_enhancers) <- headers[2:length(headers)]


# Variables de composition nucléotidiques pour les TSS
quadri_proms = read.table(paste0("Compo_nucleo_promoters/genes_", ct, "_", as.character(nb_nucl),".txt"), sep = ';', h = T, na.strings="")
headers = read.table(paste0("Compo_nucleo_promoters/genes_", ct, "_", as.character(nb_nucl), ".txt"), sep = ';', h = F, na.strings="")[1,]
row.names(quadri_proms) = quadri_proms$X
quadri_proms = quadri_proms[ ,(names(quadri_proms) != 'X')]
colnames(quadri_proms) <- headers[2:length(headers)]

#Jeu de données comprenant les paires enhancers-TSS positives et négatives
pairs = read.table(paste0("Length_AUC_plot/EG_pairs_", ct, ".txt"), sep = ';', h = T)
cols = c("pairs", paste0("E_", row.names(quadri_proms)), paste0("P_", row.names(quadri_proms)), "chr", "position", "distance", "interaction")
data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
colnames(data) <- cols
data$interaction = pairs$Interaction
data$distance = pairs$Distance
data$position = pairs$Enhancer_center
data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]

#Pour chaque paire, leur attribuer les variables de taux de nucléotides
for(i in 1:length(pairs$Enhancer)){
  data$pairs[i] = paste0(pairs$Enhancer[i],"/", pairs$Gene[i])
  data[i, paste0("E_", row.names(quadri_enhancers))] = quadri_enhancers[, pairs$Enhancer[i]]
  data[i, paste0("P_", row.names(quadri_proms))] = quadri_proms[, pairs$Gene[i]]
}
#order by chromosomal position
data <- data[order(data$pairs),]
rownames(data) <- data$pairs
data <- na.omit(data)

par(mfrow = c(3,4))
for(var in row.names(quadri_enhancers)){
  hist(data[,paste0('E_',var)], breaks = 30, col = rgb(1,0,0,0.3), ylim = c(0,100), main = var)
  hist(data[,paste0('P_',var)], breaks = 30, col = rgb(0,0,1,0.3), add = T)
}
hist(data[data$interaction==1,]$E_AT, breaks = 100, col = rgb(1,0,0,0.3))
hist(data[data$interaction==1,]$P_AT, breaks = 100, col = rgb(0,0,1,0.3), add = T)

#Separate training set and evaluation set
chrs = unique(data$chr)
test_pairs = c()
chosen_enhancers = c()
while(length(test_pairs)< round(1/5*length(data$pairs))){
  if(length(chrs > 0)){
    chr <- sample(chrs, size = 1)
    chrs = chrs[chrs != chr]
  }
  else{
    chrs = unique(data$chr)
    chr <- sample(chrs, size = 1)
  }
  data_chr <- subset(data, data$chr == chr)
  print(chr)
  #on ne prend pas un enhancer deja dans le jeu de test
  data_chr <- subset(data_chr, !str_split_fixed(data_chr$pairs, '/', 2)[,1] %in% chosen_enhancers)
  enh <- sample(unique(str_split_fixed(data_chr$pairs, '/', 2)[,1]), size = 1)
  print(enh)
  chosen_enhancers = c(chosen_enhancers, enh)
  data_chr_enh <- subset(data_chr, str_split_fixed(data_chr$pairs, '/', 2)[,1] == enh)
  print(data_chr_enh$pairs)
  test_pairs = c(test_pairs, data_chr_enh$pairs)
}
print(test_pairs)

data_test <- data[test_pairs,]
head(data_test)
table(data_test$chr)
data_train <- data[!rownames(data) %in% rownames(data_test),]

#naive way, with a potential biais?
data_train <- data[1:round(dim(data)[1]*4/5),]
data_test <- data[which(!rownames(data) %in% rownames(data_train)),]


#defining variables
variables_ep = names(data)[which(names(data)!="interaction" & names(data)!="pairs" & names(data)!="chr"
                              & names(data)!="position")]

variables_e = variables_ep[which(!variables_ep %in% paste0("E_", row.names(quadri_proms)))]


#model with enhancers and promoter variables
x_ep <- as.matrix(data_train[, variables_ep])
glm_ep <- cv.glmnet(x_ep, y, alpha = 1, family="binomial", keep = T, nfolds = 10)

#model with enhancer variables
x_e  <- as.matrix(data_train[, variables_e])
glm_e <- cv.glmnet(x_e, y, alpha = 1, family="binomial", keep = T, nfolds = 10)

#evaluate the models and compare them
x_test_ep <- as.matrix(data_test[, variables_ep])
x_test_e <- as.matrix(data_test[, variables_e])

y_test <- as.matrix(data_test$interaction)

pred_ep <- as.numeric(predict(glm_ep, x_test_ep, type = "response"))
pred_e <- as.numeric(predict(glm_e, x_test_e, type = "response"))


iLambdaMinq_ep = which(glm_ep$lambda==glm_ep$lambda.min)
iLambdaMinq_e = which(glm_e$lambda==glm_e$lambda.min)

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
auc_e <- roc(data_test$interaction, pred_e)$auc

curve_ep = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = pred_ep)
curve_e = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = pred_e)


par(mfrow = c(1,1))
png(paste0('Figures/roc',ct, "_", as.character(nb_nucl),'.png'))
plot(curve_ep[2,], curve_ep[1,], type = "l", col = "blue", 
     main = paste("ROC", ct, "AUC : ", as.character(round(auc_ep,4))), xlim = c(0,1), ylim = c(0,1),
     xlab = "FP rate", ylab = "TP rate")
lines(curve_e[2,], curve_e[1,], type = "l", col = "green")
legend(0.3, 0.2, legend=c(paste0("Enhancers and promoters features, auc = ", as.character(round(auc_ep,4))), 
                          paste0("Enhancers features, auc = ", as.character(round(auc_e,4))), "Pormoter features"),
       col=c("blue", "red", "green"), lty=1, cex=0.8)
abline(c(0,0),c(1,1))
dev.off()

varselq =names(which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0))
varselq
varsel_e
varsel_e = names(which(glm_e$glmnet.fit$beta[,iLambdaMinq_e]!=0))
glm_ep$glmnet.fit$beta[,iLambdaMinq_ep][which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0)]

table(data$chr)
table(data_test$chr)
#data_description <- function(data_train){
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
  
  print(length(which(enh_test %in% enh_train)))
  print(length(which(prom_test %in% prom_train)))
  print(prom_test[which(prom_test %in% prom_train)])
}
train_test_intersect()

data_description(data_train)
