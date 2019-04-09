rm(list=ls())
setwd("~/DataFantom")
library(gridExtra)
source("make_datasets.R")
source("make_models.R")

ct = commandArgs(trailingOnly=TRUE)[1]

pairs = read.table(paste0("EG_pairs/EG_pairs_", ct, ".txt"), sep = ';', h = T)
dexter_vars_train_enh <- read.table(paste0("Dexter_results/Regression_gene_expr/", ct,"/models/all_domains.dat_training_set.log.matrix"), h = T, sep = ' ')
dexter_vars_test_enh <- read.table(paste0("Dexter_results/Regression_gene_expr/", ct,"/models/all_domains.dat_testing_set.log.matrix"), h = T, sep = ' ')

dexter_vars_enh <- rbind.data.frame(dexter_vars_train_enh, dexter_vars_test_enh, make.row.names = TRUE )
row.names(dexter_vars_enh) = dexter_vars_enh$sequence
dexter_vars_enh = dexter_vars_enh[ ,(names(dexter_vars_enh) != 'sequence' & names(dexter_vars_enh)!= "expression")]

cols = c("pairs", paste0("E_", names(dexter_vars_enh)), "chr", "position", "distance", "interaction")
data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
colnames(data) <- cols
data$interaction = pairs$Interaction
data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]
data$gene_expr = pairs$gene_expression

for(i in 1:length(pairs$Enhancer)){
  data$pairs[i] = paste0(pairs$Enhancer[i],"/", pairs$Gene[i])
  data[i, paste0("E_", names(dexter_vars_enh))] = dexter_vars_enh[data$pairs[i],]
}

x  <- as.matrix(data[, paste0("E_", names(dexter_vars_enh))])
y <- as.matrix(data$gene_expr)
glm <- cv.glmnet(x, y, alpha = 1, family="gaussian", keep = T, nfolds = 10)


pred <- as.numeric(predict(glm, x))
error <- abs(y-pred)
data$error=error

test <- ks.test(as.vector(data[data$interaction==1,]$error), as.vector(data[data$interaction==0,]$error))
test$p.value

png(paste0('Figures/error_hist_pos_neg_',ct,'.png'))
hist(data[data$interaction==1,]$error, col = rgb(0,1,0,0.3), breaks = 25,
     main = paste("Distribution de l'ereur en regression, pval du ks : ", test$p.value),
     xlab = "error")
hist(data[data$interaction==0,]$error, col = rgb(1,0,0,0.3), breaks = 25, add=T)
dev.off()

summary(data[data$interaction==1,]$error)
summary(data[data$interaction==0,]$error)