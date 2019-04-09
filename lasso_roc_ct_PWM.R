library(glmnet)
library(stringr)
rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)

ct = commandArgs(trailingOnly=TRUE)[1]
nb_nucl = commandArgs(trailingOnly=TRUE)[2]
load("Data/unique_annot_ct.RData")


overlapping_enh <- function(enh, others){
  start = as.numeric(str_split_fixed(str_split_fixed(enh, '-', 2)[,1], ':', 2)[,2])
  end = as.numeric(str_split_fixed(enh, '-', 2)[,2])
  starts = as.numeric(str_split_fixed(str_split_fixed(others, '-', 2)[,1], ':', 2)[,2])
  ends = as.numeric(str_split_fixed(others, '-', 2)[,2])
  overlapping_end = others[starts <= end & starts >= start]
  overlapping_start = others[ends >= start & ends <= end]
  return(union(overlapping_start, overlapping_end))
}


# Variables de composition nucléotidiques pour les enhancers

ct = "CNhs11790"
nb_nucl = 2
# Variables de composition nucléotidiques pour les enhancers
quadri_enhancers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl),".txt"), sep = ';', h = T, na.strings="")
headers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl), ".txt"), sep = ';', h = F, na.strings="")[1,]
row.names(quadri_enhancers) = quadri_enhancers$X
quadri_enhancers = quadri_enhancers[ ,(names(quadri_enhancers) != 'X')]
colnames(quadri_enhancers) <- headers[2:length(headers)]

pairs = read.table(paste0("EG_pairs/EG_pairs_", ct, ".txt"), sep = ';', h = T)

#pwm_enh <- read.table("MotifSearch_enhancers/mat_score_enhancers.txt", h = T, sep = ',', row.names = 1, check.names = FALSE)
#save(headers, file = "Data/pwm_enhancers_headers.RData")
#save(pwm_enh, file = "Data/pwm_enhancers.RData")

load(file = "Data/pwm_enhancers.RData")

pwm_ct <- pwm[,names(pwm) %in% pairs$Enhancer]

rownames(pwm)
#Jeu de données comprenant les paires enhancers-TSS positives et négatives

cols = c("pairs", paste0("E_", row.names(quadri_enhancers)), paste0("E_", rownames(pwm_ct)),"chr", "position", "distance", "interaction")
data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
colnames(data) <- cols
data$interaction = pairs$Interaction
data$distance = pairs$Distance
data$position = pairs$Enhancer_center
data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]

#Pour chaque paire, leur attribuer les variables de taux de nucléotides
for(i in 1:length(pairs$Enhancer)){
  data$pairs[i] = paste0(pairs$Enhancer[i],"/", pairs$Gene[i])
  data[i, paste0("E_", rownames(pwm_ct))] = pwm_ct[,pairs$Enhancer[i]]
  data[i, paste0("E_", row.names(quadri_enhancers))] = quadri_enhancers[, pairs$Enhancer[i]]
}
#order by chromosomal position
data <- data[order(data$pairs),]
rownames(data) <- data$pairs
data <- na.omit(data)

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
  data_chr <- data[data$chr == chr, ]
  #on ne prend pas un enhancer deja dans le jeu de test
  data_chr <- data_chr[!str_split_fixed(data_chr$pairs, '/', 2)[,1] %in% chosen_enhancers,]
  #print(unique(str_split_fixed(data_chr$pairs, '/', 2)[,1]))
  if(length(unique(str_split_fixed(data_chr$pairs, '/', 2)[,1]))>1){
    enh <- sample(unique(str_split_fixed(data_chr$pairs, '/', 2)[,1]), size = 1)
    overlapping_enhs <- overlapping_enh(enh, str_split_fixed(data_chr$pairs, '/', 2)[,1])
    chosen_enhancers = c(chosen_enhancers, overlapping_enhs)
    data_chr_enh <- data_chr[str_split_fixed(data_chr$pairs, '/', 2)[,1] %in% overlapping_enhs,]
    test_pairs = c(test_pairs, data_chr_enh$pairs)
  }
  
}
data_test <- data[test_pairs,]
data_train <- data[!rownames(data) %in% rownames(data_test),]



#defining variables
variables_ep = names(data)[which(names(data)!="interaction" & names(data)!="pairs" & names(data)!="chr"
                                 & names(data)!="position")]

variables_e = variables_ep[which(!variables_ep %in% paste0("E_", rownames(pwm_ct)))]

y <- as.matrix(data_train$interaction)

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
library(ggplot2)

png(paste0('Figures/roc',ct, '_pwm_dinucl.png'))
roc <- data.frame(x = curve_ep[2,], y = curve_ep[1,])
a <- ggplot(data = roc, aes(x = x, y = y))
a + geom_line(aes(colour = paste0("PWM + nucleotides AUC : ", round(auc_ep,3)))) + 
  geom_line(aes(x = seq(0,1,by = 1/(dim(roc)[1]-1)), y = seq(0,1,by =  1/(dim(roc)[1]-1)),colour = "Random")) + 
  geom_line(aes(x = curve_e[2,], y = curve_e[1,],colour = paste0("nucleotidesAUC : ", round(auc_e,3)))) + 
  xlab("FP rate") + ylab("TP rate")  + theme(legend.position = c(0.6,0.3), legend.justification = c(0,1)) + 
  scale_color_manual(name = "Legend",values = c( "purple", "green","black"))+ 
  ggtitle(paste("ROC", ct, ' (', unique_sample_annot[ct,]$lib_id, ')')) 
dev.off()


png(paste0('Figures/roc',ct, '_dexter.png'))
plot(curve_ep[2,], curve_ep[1,], type = "l", col = "blue", 
     main = paste("ROC", ct, " PWM AUC : ", as.character(round(auc_ep,4))), xlim = c(0,1), ylim = c(0,1),
     xlab = "FP rate", ylab = "TP rate")
#lines(curve_e[2,], curve_e[1,], type = "l", col = "red")
#legend(0.3, 0.2, legend=c(paste0("Enhancers and promoters features, auc = ", as.character(round(auc_ep,4))), 
#paste0("Enhancers features, auc = ", as.character(round(auc_e,4)))),
#col=c("blue", "red"), lty=1:2, cex=0.8)
abline(c(0,0),c(1,1))
dev.off()


varselq =names(which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0))
glm_ep$glmnet.fit$beta[,iLambdaMinq_ep][which(glm_ep$glmnet.fit$beta[,iLambdaMinq_ep]!=0)]


#data_description(data)
