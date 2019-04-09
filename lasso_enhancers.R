library(glmnet)
library(stringr)
library(pROC) 

rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)

#essayer de prédire que les enhancers

# Variables de composition nucléotidiques pour les enhancers
quadri_enhancers = read.table("quadri_enhancers_CNhs11745.txt", sep = ';', h = T, na.strings="")
headers = read.table("quadri_enhancers_CNhs11745.txt", sep = ';', h = F, na.strings="")[1,]
row.names(quadri_enhancers) = quadri_enhancers$X
quadri_enhancers = quadri_enhancers[ ,(names(quadri_enhancers) != 'X')]
colnames(quadri_enhancers) <- headers[2:length(headers)]

#Jeu de données comprenant les paires enhancers-TSS positives et négatives
pairs = read.table("EG_pairs_CNhs11745.txt", sep = ';', h = T)
head(pairs)

cols = c("pairs", paste0("E_", row.names(quadri_enhancers)), "chr", "position", "distance", "interaction")
data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
colnames(data) <- cols

data$interaction = pairs$Interaction
data$distance = pairs$EG.distance
data$position = pairs$Enhancer_center
data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]

#Pour chaque paire, leur attribuer les variables de taux de nucléotides
for(i in 1:length(pairs$Enhancer)){
  data$pairs[i] = paste0(pairs$Enhancer[i],"-", pairs$Gene[i])
  data[i, paste0("E_", row.names(quadri_enhancers))] = quadri_enhancers[, pairs$Enhancer[i]]
}

rownames(data) <- data$pairs
data <- na.omit(data)

summary(pairs$Enhancer)
summary(pairs)

#Séparer les données en jeu de test et d'apprentissage après tri par ordre chromosomique
# (dataset dejà trié)

pos = subset(data, data$interaction==1)
neg = subset(data, data$interaction==0)
data_train <- rbind.data.frame(pos[1:round(dim(pos)[1]*4/5),], neg[1:round(dim(neg)[1]*4/5),])
data_test <- data[which(!rownames(data) %in% rownames(data_train)),]

variables = names(data)[which(names(data)!="interaction" & names(data)!="pairs" & names(data)!="chr"
                              & names(data)!="position")]
x <- as.matrix(data_train[, variables])
y <- as.matrix(data_train$interaction)
glmq <- cv.glmnet(x, y, alpha = 1, family="binomial", keep = T, nfolds = 10)


iLambdaMinq = which(glmq$lambda==glmq$lambda.min)
#predq = glmq$fit.preval[,iLambdaMinq]


#predictions uniquement sur le jeu de test, pour éviter les biais
x_test <- as.matrix(data_test[, variables])
y_test <- as.matrix(data_test$interaction)
predq <- as.numeric(predict(glmq, x_test, type = "response"))

positivesrates = function(seuil, d, pred){
  pbin = ifelse(pred>seuil, 1, 0)
  err = d$interaction - pbin
  err2 = d$interaction + pbin
  TP = sum(err2==2)
  FP = sum(err == -1)
  c(TP/sum(d$interaction==1), FP/sum(d$interaction==0))
}
par(mfrow = c(1,1))
a = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = predq)

plot(a[2,], a[1,], type = "l", col = "blue", main = paste("ROC enhancers", "CNhs11745"), xlim = c(0,1), ylim = c(0,1))
abline(c(0,0),c(1,1))



varselq =names(which(glmq$glmnet.fit$beta[,iLambdaMinq]!=0))
varselq
 

roc(data_test$interaction, predq)$auc

