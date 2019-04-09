library(glmnet)
library(stringr)
library(pROC) 

rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)


nb_nucl = 2
cts <- read.table("cell_types_plt.txt", h = F, sep = '\t')
cts = cts[1:length(cts)-1]

overlapping_enh <- function(enh, others){
  start = as.numeric(str_split_fixed(str_split_fixed(enh, '-', 2)[,1], ':', 2)[,2])
  end = as.numeric(str_split_fixed(enh, '-', 2)[,2])
  starts = as.numeric(str_split_fixed(str_split_fixed(others, '-', 2)[,1], ':', 2)[,2])
  ends = as.numeric(str_split_fixed(others, '-', 2)[,2])
  overlapping_end = others[starts <= end & starts >= start]
  overlapping_start = others[ends >= start & ends <= end]
  return(union(overlapping_start, overlapping_end))
}

positivesrates = function(seuil, d, pred){
  pbin = ifelse(pred>seuil, 1, 0)
  err = d$interaction - pbin
  err2 = d$interaction + pbin
  TP = sum(err2==2)
  FP = sum(err == -1)
  c(TP/sum(d$interaction==1), FP/sum(d$interaction==0))
}


get_auc_dexter <- function(ct){
  dexter_vars_train <- read.table(paste0("Dexter_results/", ct,"/models/all_domains.dat_training_set.log.matrix"), h = T, sep = ' ')
  dexter_vars_test <- read.table(paste0("Dexter_results/", ct,"/models/all_domains.dat_testing_set.log.matrix"), h = T, sep = ' ')
  
  dexter_vars <- rbind.data.frame(dexter_vars_train, dexter_vars_test, make.row.names = TRUE )
  
  row.names(dexter_vars) = dexter_vars$sequence
  dexter_vars = dexter_vars[ ,(names(dexter_vars) != 'sequence' & names(dexter_vars)!= "expression")]
  
  #Jeu de données comprenant les paires enhancers-TSS positives et négatives
  pairs = read.table(paste0("Length_AUC_plot/EG_pairs_", ct, ".txt"), sep = ';', h = T)
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
  
  #variables_e = variables_ep[which(!variables_ep %in% paste0("P_", row.names(quadri_proms)))]
  
  y <- as.matrix(data_train$interaction)
  
  #model with enhancers variables
  x_ep <- as.matrix(data_train[, variables_ep])
  glm_ep <- cv.glmnet(x_ep, y, alpha = 1, family="binomial", keep = T, nfolds = 10)
  

  #evaluate the models and compare them
  x_test_ep <- as.matrix(data_test[, variables_ep])
  y_test <- as.matrix(data_test$interaction)
  pred_ep <- as.numeric(predict(glm_ep, x_test_ep, type = "response"))
  return(roc(data_test$interaction, pred_ep)$auc)
}

get_auc_dinucl <- function(ct){
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
  cols = c("pairs", paste0("E_", row.names(quadri_enhancers)), paste0("P_", row.names(quadri_enhancers)), "chr", "position", "distance", "interaction")
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
  
  variables_e = variables_ep[which(!variables_ep %in% paste0("P_", row.names(quadri_proms)))]
  y <- as.matrix(data_train$interaction)
  
  #model with enhancer variables
  x_e  <- as.matrix(data_train[, variables_e])
  glm_e <- cv.glmnet(x_e, y, alpha = 1, family="binomial", keep = T, nfolds = 10)
  
  #evaluate the models and compare them
  x_test_e <- as.matrix(data_test[, variables_e])
  y_test <- as.matrix(data_test$interaction)
  
  pred_e <- as.numeric(predict(glm_e, x_test_e, type = "response"))
  return(c(dim(data)[1], roc(data_test$interaction, pred_e)$auc))
}

dexter = c()
dinucl = c()
for(ct in cts){
  print(ct)
  res = get_auc_dinucl(ct)
  dinucl = c(dinucl, res[2])
  dexter = c(dexter, get_auc_dexter(ct))
}


png(paste0('Figures/dexter_dinucl_100.png'))
plot(x = dinucl, y = dexter, pch = 16, col = 'darkblue',
     main = 'model AUCs of dexter variables against dinucleotide variables')
abline(c(0,0),c(1,1))
dev.off()
