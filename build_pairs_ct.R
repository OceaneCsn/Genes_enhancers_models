library(stringr)
options(stringsAsFactors=F)
setwd("~/DataFantom")

#Expression matrices
load("Data/fantom.prom.tpm.RData")
load("Data/fantom.enh.tpm.filtered.RData")
load("Data/fantom.sample.annot.RData")


Me_fpkm <- Me
cell_types = colnames(Me_fpkm)
ct = commandArgs(trailingOnly=TRUE)[1]

mode = commandArgs(trailingOnly=TRUE)[2]

#ct = "CNhs11889"
#fantom associations based on expression correlation
associations <- read.table("Data/hg19_enhancer_promoter_correlations_distances_cell_type.txt", sep = '\t', h = T)

#fantom annotations to choose only TSS associated with coding genes
annot <- read.table("Data/annotation_tss_peaks.txt", h = T, sep = '\t')
annot <- na.omit(annot)
annot_coding <- annot[which(str_split_fixed(annot$short_description, '@', 2)[,1]=="p1"),]
Mg_coding <- Mg_fpkm[rownames(Mg_fpkm) %in% annot_coding$X00Annotation,]

#log expression to compare and have better distributions
Mg_coding_log <- log(Mg_coding+min(Mg_coding[Mg_coding!=0]))

#choose one sample per cell type
unique_sample_annot <- sample.annot[!duplicated(sample.annot$lib_id),]
unique_sample_annot$lib_id <- as.character(unique_sample_annot$lib_id)



# Computes the number of EG pairs found in the 600 distinct cell types of the fantom dataset
stats_pairs_ct <- function(){
  nb_pairs = data.frame(row.names = rownames(unique_sample_annot))
  nb_pairs$pairs = rep(0, length(rownames(nb_pairs)))
  for(cell_t in rownames(unique_sample_annot)){
    expressed_proms <- names(Mg_coding[,cell_t][which(Mg_coding[,cell_t]>0)])
    expressed_enh <- names(Me_fpkm[,cell_t][which(Me_fpkm[,cell_t]>0)])
    print(length(expressed_enh))
    pos = associations[associations$enhancer %in% expressed_enh,]
    pos = pos[pos$promoter %in% expressed_proms,]
    nb_pairs[cell_t, "pairs"] = dim(pos)[1]
  }
  return(nb_pairs)
}
#par(mfrow = c(1,1))
#counts_ct <- stats_pairs_ct()
#hist(counts_ct$pairs, breaks = 100)
#counts_ct <- counts_ct[order(counts_ct$pairs, decreasing = T),]
#counts_ct$description <- unique_sample_annot[rownames(counts_ct), "lib_id"]


#builds a balanced dataset for the chosen cell type, giving E-G pairs, their distance,
#and weather their interact or not
build_dataset <- function (ct, folder, gene_expr = T){
  
  expressed_proms <- names(Mg_coding[,ct][which(Mg_coding[,ct]>0)])
  expressed_enh <- names(Me_fpkm[,ct][which(Me_fpkm[,ct]>0)])
  
  ################  positive sequences active enhancers and active TSS for ct : 
  
  pos = associations[associations$enhancer %in% expressed_enh,]
  pos = pos[pos$promoter %in% expressed_proms,c("enhancer", "promoter", "distance")]
  pos$Interaction = 1
  colnames(pos) <-c("Enhancer", "Gene", "Distance", "Interaction")
  
  ################## negative sequences : unexpressed enhancers with active TSS 
  
  unexpressed_enh <- names(Me_fpkm[,ct][which(Me_fpkm[,ct]==0)])
  unexpressed_enh <- sample(unexpressed_enh, size = length(unexpressed_enh), replace = F)
  
  neg = associations[associations$enhancer %in% unexpressed_enh,]
  neg = neg[neg$promoter %in% expressed_proms, c("enhancer", "promoter", "distance")]
  
  # control the promoter expression so the distribution of expression is the same for the
  # promoters in the negative pairs and the positive pairs
  controled_proms = c()
  
  #if there are more negatives than positives
  if(dim(neg)[1]> dim(pos)[1]){
    posproms = pos$Gene
    negproms = neg[neg$promoter %in% expressed_proms,]$promoter
    negproms <- sample(negproms, size = length(negproms), replace = F)
    while(length(controled_proms) < 4/5*length(pos$Gene)){
      posprom = posproms[1]
      seuil = 0
      while(length(which(abs(Mg_coding_log[as.character(negproms),ct]-Mg_coding_log[as.character(posprom),ct]) < seuil))==0){
        if(seuil>5){
          seuil = seuil+0.5
        }
        else{
          seuil = seuil +0.05
        }
      }
      negprom = names(which(abs(Mg_coding_log[as.character(negproms),ct]-Mg_coding_log[as.character(posprom),ct]) < seuil))[1]
      index = match(negprom, negproms)
      negproms <- negproms[-index]
      controled_proms = c(controled_proms, negprom)
      if(length(posproms)>1){
        posproms <- posproms[2:length(posproms)]
      }
      else{posproms = c()}
    }
    neg = neg[neg$promoter %in% controled_proms, c("enhancer", "promoter", "distance")]
    neg = neg[sample(rownames(neg), size = min(length(pos$Enhancer), length(neg$promoter))),]
    neg$Interaction = 0
    colnames(neg) <-c("Enhancer", "Gene", "Distance", "Interaction")
    
  }
  #if there are more positives than negatives
  else{
    controled_enhs = c()
    negproms = neg$promoter
    pos = pos[sample(rownames(pos), size = length(rownames(pos))),]
    posproms = pos$Gene
    posenhs = pos$Enhancer
    while(length(controled_proms) < length(neg$promoter)){
      negprom = negproms[1]
      seuil = 0
      while(length(which(abs(Mg_coding_log[as.character(posproms),ct]-Mg_coding_log[as.character(negprom),ct]) < seuil))==0){
        if(seuil>5){
          seuil = seuil+0.5
        }
        else{
          seuil = seuil +0.05
        }
      }
      posprom = names(which(abs(Mg_coding_log[as.character(posproms),ct]-Mg_coding_log[as.character(negprom),ct]) < seuil))[1]
      index = match(posprom, posproms)
      posenh = posenhs[index]
      posproms <- posproms[-index]
      posenhs <- posenhs[-index]
      controled_proms = c(controled_proms, posprom)
      controled_enhs = c(controled_enhs, posenh)
      if(length(negproms)>1){
        negproms <- negproms[2:length(negproms)]
      }
      else{negproms = c()}
    }
    pos = pos[pos$Gene %in% controled_proms, c("Enhancer", "Gene", "Distance", "Interaction")]
    pos =  pos[pos$Enhancer %in% controled_enhs, c("Enhancer", "Gene", "Distance", "Interaction")]
    neg$interaction = 0
    colnames(neg) <-c("Enhancer", "Gene", "Distance", "Interaction")
    neg <- neg[c("Enhancer", "Gene", "Distance", "Interaction")]
  }
  par(mfrow = c(1,2))
  hist(Mg_coding_log[as.character(pos$Gene),ct], breaks = 1000)
  hist(Mg_coding_log[as.character(neg$Gene),ct], breaks = 1000)
  par(mfrow = c(1,1))
  data_train <- rbind.data.frame(pos, neg, make.row.names = TRUE )
  if(gene_expr){
    data_train$gene_expression = 0
    for(i in rownames(data_train)){
      data_train[i, "gene_expression"] = Mg_coding_log[data_train[i, "Gene"], ct]
    }
  }
  
  write.table(data_train, file = paste0(folder, "/EG_pairs_",ct,".txt"), row.names=FALSE, quote = F, sep = ';')
  return(data_train)
}

stats_dataset_00 <- function (){
  
  cts <- sample(rownames(unique_sample_annot), size = 500)
  
  sizes_00 <- c()
  sizes_11 <- c()
  sizes_01 <- c()
  sizes_10 <- c()
  
  for (ct in cts){
    expressed_proms <- names(Mg_coding[,ct][which(Mg_coding[,ct]>0)])
    expressed_enh <- names(Me_fpkm[,ct][which(Me_fpkm[,ct]>0)])
    
    ################  positive sequences 
    
    pos_11 = associations[associations$promoter %in% expressed_proms & associations$enhancer %in% expressed_enh, ]
    
    pos_00 = associations[!associations$promoter %in% expressed_proms & !associations$enhancer %in% expressed_enh,]
    pos_00 = pos_00[pos_00$promoter %in% rownames(Mg_coding),]
    
    ################## negative sequences      
    
    neg_01 = associations[!associations$enhancer %in% expressed_enh & associations$promoter %in% expressed_proms,]
    
    neg_10 = associations[associations$enhancer %in% expressed_enh & !associations$promoter %in% expressed_proms,]
    neg_10 = neg_10[neg_10$promoter %in% rownames(Mg_coding),]
    
    sizes_00 <- c(sizes_00, dim(pos_00)[1])
    sizes_11 <- c(sizes_11, dim(pos_11)[1])
    sizes_01 <- c(sizes_01, dim(neg_01)[1])
    sizes_10 <- c(sizes_10, dim(neg_10)[1])
  }
  
  sizes <- c(sizes_00, sizes_01, sizes_11, sizes_10)
  types <- c(rep('00', length(sizes_00)), rep('01', length(sizes_00)), rep('11', length(sizes_00)), rep('10', length(sizes_00)))
  
  stats <- data.frame('sizes' = sizes, 'types' = types)
  
  
  stats_ct <- data.frame('00' = sizes_00, '01' = sizes_01, '11' = sizes_11, '10' = sizes_10)
  rownames(stats_ct) <- cts
  
  stats_ct_filtered <- stats_ct[stats_ct$X10 < stats_ct$X11 & stats_ct$X10 < stats_ct$X00 & stats_ct$X10 < stats_ct$X01, ]
  stats_ct_filtered <- stats_ct_filtered[stats_ct$X10 > 250,]
  
  stats_ct_filtered <- stats_ct_filtered[stats_ct_filtered$X11 < stats_ct_filtered$X01,]
  
  ggplot(stats, aes(factor(types),sizes)) +
    geom_violin(aes(fill = factor(types)))
  
}



build_dataset_00 <- function (ct){
  
  
  expressed_proms <- names(Mg_coding[,ct][which(Mg_coding[,ct]>0)])
  expressed_enh <- names(Me_fpkm[,ct][which(Me_fpkm[,ct]>0)])
  
  ################  positive sequences 
  
  pos_11 = associations[associations$promoter %in% expressed_proms & associations$enhancer %in% expressed_enh, ]
  
  pos_00 = associations[!associations$promoter %in% expressed_proms & !associations$enhancer %in% expressed_enh,]
  pos_00 = pos_00[pos_00$promoter %in% rownames(Mg_coding),]
  
  ################## negative sequences      
  
  neg_01 = associations[!associations$enhancer %in% expressed_enh & associations$promoter %in% expressed_proms,]
  
  neg_10 = associations[associations$enhancer %in% expressed_enh & !associations$promoter %in% expressed_proms,]
  neg_10 = neg_10[neg_10$promoter %in% rownames(Mg_coding),]
  
  
  ############### sampling 00,11 and 01 so they have the same size as 10
  
  pos_00 = pos_00[sample(rownames(pos_00), size = dim(neg_10)[1]),]
  
  # for 11 and 01, choosing promoters with the same expression
  
  pos_11_ech = data.frame(matrix(ncol = length(colnames(pos_11)), nrow = dim(pos_11)[1]))
  colnames(pos_11_ech) <- colnames(pos_11)
  
  neg_01_ech = data.frame(matrix(ncol = length(colnames(pos_11)), nrow = dim(pos_11)[1]))
  colnames(neg_01_ech) <- colnames(pos_11)
  
  pos_11_temp = data.frame(pos_11)
  
  rownames(neg_01) <- seq(1:dim(neg_01)[1])
  neg_01_temp = data.frame(neg_01)
  
  for (row in rownames(pos_11_ech)){
    
    pos_row = pos_11_temp[sample(rownames(pos_11_temp), size = 1),]
    posprom = pos_row$promoter
    
    pos_11_temp = pos_11_temp[rownames(pos_11_temp) != rownames(pos_row),]
    pos_11_ech[row,] = pos_row
    
    seuil = 0
    while(length(which(abs(Mg_coding_log[as.character(neg_01$promoter),ct]-Mg_coding_log[as.character(posprom),ct]) < seuil))==0){
      if(seuil>1){
        seuil = seuil+0.05
      }
      else{
        seuil = seuil+0.01
      }
    }
    neg_row = neg_01_temp[as.vector(which(abs(Mg_coding_log[as.character(neg_01_temp$promoter),ct]-Mg_coding_log[as.character(posprom),ct]) < seuil))[1],]
    neg_01_temp = neg_01_temp[rownames(neg_01_temp) != rownames(neg_row),]
    neg_01_ech[row,] = neg_row
  }
  
  pos_11_ech$pair_type = '11'
  pos_00$pair_type = '00'
  neg_01_ech$pair_type = '01'
  neg_10$pair_type = '10'
  pos <- rbind.data.frame(pos_11_ech, pos_00)
  pos$interaction = 1
  neg <- rbind.data.frame(neg_10, neg_01_ech)  
  neg$interaction = 0
  
  data <- rbind.data.frame(pos, neg)
  
  data <- data[c('enhancer', 'promoter', 'distance', 'interaction','pair_type')]
  colnames(data) <- c("Enhancer", "Gene", "Distance", "Interaction", "Pair_type")
  
  data <- na.omit(data)
  
  data$gene_expression = 0
  for(i in rownames(data)){
    print(Mg_coding_log[data[i, "Gene"], ct])
    data[i, "gene_expression"] = Mg_coding_log[data[i, "Gene"], ct]
  }
  
  write.table(data, file = paste0("EG_pairs/EG_pairs_",ct,".txt"), row.names=FALSE, quote = F, sep = ';')
  
  return(data)
}

if(mode=='00'){
  print('1100')
  data_train <- build_dataset_00(ct)
}
if(mode=='11'){
  print('11')
  data_train <- build_dataset(ct, 'EG_pairs')
}