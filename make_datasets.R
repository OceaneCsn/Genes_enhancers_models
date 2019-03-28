library(glmnet)
library(stringr)
setwd("~/DataFantom")
options(stringsAsFactors=F)

create_dataset <- function(ct, variables = c("nucl")){
  nb_nucl = 2
  vars = c()
  
  #Jeu de données comprenant les paires enhancers-TSS positives et négatives
  pairs = read.table(paste0("EG_pairs/EG_pairs_", ct, ".txt"), sep = ';', h = T)
  
  if("nucl" %in% variables){
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
    vars <- c(paste0("E_", row.names(quadri_proms)), paste0("P_", row.names(quadri_proms)))
  }
  
  if("dexter" %in% variables){
    #Variables from dexter exploration for enhancers
    dexter_vars_train_enh <- read.table(paste0("Dexter_results/Enhancers/", ct,"/models/all_domains.dat_training_set.log.matrix"), h = T, sep = ' ')
    dexter_vars_test_enh <- read.table(paste0("Dexter_results/Enhancers/", ct,"/models/all_domains.dat_testing_set.log.matrix"), h = T, sep = ' ')
    dexter_vars_enh <- rbind.data.frame(dexter_vars_train_enh, dexter_vars_test_enh, make.row.names = TRUE )
    
    row.names(dexter_vars_enh) = dexter_vars_enh$sequence
    dexter_vars_enh = dexter_vars_enh[ ,(names(dexter_vars_enh) != 'sequence' & names(dexter_vars_enh)!= "expression")]
    vars = c(vars, paste0('E_',names(dexter_vars_enh)))
    
    #Variables from dexter exploration for enhancers
    dexter_vars_train_prom <- read.table(paste0("Dexter_results/Promoters/", ct,"/models/all_domains.dat_training_set.log.matrix"), h = T, sep = ' ')
    dexter_vars_test_prom <- read.table(paste0("Dexter_results/Promoters/", ct,"/models/all_domains.dat_testing_set.log.matrix"), h = T, sep = ' ')
    dexter_vars_prom <- rbind.data.frame(dexter_vars_train_prom, dexter_vars_test_prom, make.row.names = TRUE )
    
    row.names(dexter_vars_prom) = dexter_vars_prom$sequence
    dexter_vars_prom = dexter_vars_prom[ ,(names(dexter_vars_prom) != 'sequence' & names(dexter_vars_prom)!= "expression")]
    vars = c(vars, paste0('P_',names(dexter_vars_prom)))
  }
  
  if("pwm" %in% variables){
    load(file = "Data/pwm_enhancers.RData")
    pwm_ct_enh <- pwm_enh[,names(pwm_enh) %in% pairs$Enhancer]
    vars = c(vars, paste0("E_", rownames(pwm_ct_enh)))
    
    load(file = "Data/pwm_promoters.RData")
    pwm_ct_proms <- pwm_proms[,names(pwm_proms) %in% pairs$Gene]
    vars = c(vars, paste0("P_", rownames(pwm_ct_proms)))
  }
  
  #construction de la matrice
  cols = c("pairs", vars, "chr", "position", "distance", "interaction")
  data = data.frame(matrix(ncol = length(cols), nrow = dim(pairs)[1]))
  colnames(data) <- cols
  data$interaction = pairs$Interaction
  data$distance = pairs$Distance
  data$position = pairs$Enhancer_center
  data$chr = str_split_fixed(pairs$Enhancer, ':', 2)[,1]
  
  for(i in 1:length(pairs$Enhancer)){
    data$pairs[i] = paste0(pairs$Enhancer[i],"/", pairs$Gene[i])
    if("nucl" %in% variables){
      data[i, paste0("E_", row.names(quadri_enhancers))] = quadri_enhancers[, pairs$Enhancer[i]]
      data[i, paste0("P_", row.names(quadri_proms))] = quadri_proms[, pairs$Gene[i]]
    }
    if("pwm" %in% variables){
      data[i, paste0("E_", rownames(pwm_ct_enh))] = pwm_ct_enh[,pairs$Enhancer[i]]
      data[i, paste0("P_", rownames(pwm_ct_proms))] = pwm_ct_proms[,pairs$Gene[i]]
    }
    if("dexter" %in% variables){
      #print(dexter_vars_enh[pairs$Enhancer[i],])
      data[i, paste0("E_", names(dexter_vars_enh))] = dexter_vars_enh[pairs$Enhancer[i],]
      data[i, paste0("P_", names(dexter_vars_prom))] = dexter_vars_prom[pairs$Gene[i],]
    }
  }
  print(dim(data))
  data <- data[order(data$pairs),]
  rownames(data) <- data$pairs
  data <- na.omit(data)
  return(data)
}

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

train_test_split_order <- function(data){
  #Separate training set and evaluation set
  data <- data[order(data$pairs),]
  
  data_train <- data[1:round(dim(data)[1]*4/5),]
  data_test <- data[which(!rownames(data) %in% rownames(data_train)),]
  
  return(list(data_train, data_test))
}

train_test_split <- function(data){
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
  return(list(data_train, data_test))
}

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
  par(mfrow = c(1,1))
}

train_test_intersect <- function(data_train, data_test){
  enh_train <- str_split_fixed(data_train$pairs, '/', 2)[,1]
  prom_train <- str_split_fixed(data_train$pairs, '/', 2)[,2]
  
  enh_test <- str_split_fixed(data_test$pairs, '/', 2)[,1]
  prom_test <- str_split_fixed(data_test$pairs, '/', 2)[,2]
  print('Enhencers en commun entre train et test')
  print(length(which(enh_test %in% enh_train)))
  print('Promoters en commun entre train et test')
  print(length(which(prom_test %in% prom_train)))
}