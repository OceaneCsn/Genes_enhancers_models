library(glmnet)
library(stringr)
library(pROC) 

rm(list=ls())
setwd("~/DataFantom")
options(stringsAsFactors=F)

cts <- read.table("cell_types_plt.txt", h = F, sep = '\t')
cts <- cts[-1]

nb_nucl = 2
ct = cts[10]$V1

#useful functions
overlapping_enh <- function(enh, others, already_overlapping, result){
  #recursive function to return all the overlapping enhancers of enh, ang all their overlapping ones
  start = as.numeric(str_split_fixed(str_split_fixed(enh, '-', 2)[,1], ':', 2)[,2])
  end = as.numeric(str_split_fixed(enh, '-', 2)[,2])
  starts = as.numeric(str_split_fixed(str_split_fixed(others, '-', 2)[,1], ':', 2)[,2])
  ends = as.numeric(str_split_fixed(others, '-', 2)[,2])
  res = c()
  overlapping_end = others[starts <= end & starts >= start]
  overlapping_start = others[ends >= start & ends <= end]
  already_overlapping = c(already_overlapping, enh)
  to_search = union(overlapping_start, overlapping_end)
  to_search = to_search[!to_search %in% already_overlapping]
  if(length(to_search)==0){
    return(already_overlapping)
  }
  else{
    for(en in to_search){
      result = c(result, overlapping_enh(en, others = others, already_overlapping = already_overlapping, result = result))
    }
  }
  return(unique(result))
}

positivesrates = function(seuil, d, pred){
  pbin = ifelse(pred>seuil, 1, 0)
  err = d$interaction - pbin
  err2 = d$interaction + pbin
  TP = sum(err2==2)
  FP = sum(err == -1)
  c(TP/sum(d$interaction==1), FP/sum(d$interaction==0))
}

get_auc <- function(ct){
  # Variables de composition nucléotidiques pour les enhancers
  quadri_enhancers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl),".txt"), sep = ';', h = T, na.strings="")
  headers = read.table(paste0("Compo_nucleo_enhancers/enhancers_", ct, "_", as.character(nb_nucl), ".txt"), sep = ';', h = F, na.strings="")[1,]
  row.names(quadri_enhancers) = quadri_enhancers$X
  quadri_enhancers = quadri_enhancers[ ,(names(quadri_enhancers) != 'X')]
  colnames(quadri_enhancers) <- headers[2:length(headers)]
  
  
  #Jeu de données comprenant les paires enhancers-TSS positives et négatives
  pairs = read.table(paste0("Length_AUC_plot/EG_pairs_", ct, ".txt"), sep = ';', h = T)
  cols = c("pairs", paste0("E_", row.names(quadri_enhancers)), "chr", "position", "distance", "interaction")
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
  }
  #order by chromosomal position
  data <- data[order(data$pairs),]
  rownames(data) <- data$pairs
  data <- na.omit(data)
  
  table(data$chr)
  
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
      overlapping_enhs <- overlapping_enh(enh, others = str_split_fixed(data_chr$pairs, '/', 2)[,1], already_overlapping = c(),
                                          result = c())
      chosen_enhancers = c(chosen_enhancers, overlapping_enhs)
      data_chr_enh <- data_chr[str_split_fixed(data_chr$pairs, '/', 2)[,1] %in% overlapping_enhs,]
      test_pairs = c(test_pairs, data_chr_enh$pairs)
    }
    
  }
  data_test <- data[test_pairs,]
  data_train <- data[!rownames(data) %in% rownames(data_test),]
  
  #defining variables
  variables_e = names(data)[which(names(data)!="interaction" & names(data)!="pairs" & names(data)!="chr"
                                   & names(data)!="position" & names(data)!="distance")]
  
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

aucs = c()
sizes = c()
for(ct in cts){
  print(ct)
  res = get_auc(ct)
  aucs = c(aucs, res[2])
  sizes = c(sizes, res[1])
}

# scores <- data.frame(matrix(ncol = 2, nrow = length(cts)))
# colnames(scores) <- c("AUC", "size")
# rownames(scores) = cts[1,]
# scores$AUC = aucs
# scores$size = sizes
# scores <- scores[order(scores$AUC),]



par(mfrow = c(1,1))
png(paste0('Figures/Length_AUC_110_without_distance_all_chr_in_test_true.png'))
plot(x = sizes, y = aucs, pch = 16, col = 'darkblue', ylim = c(0.53,0.86),
     main = 'AUC depending on the dataset size (number of pairs)')
dev.off()

