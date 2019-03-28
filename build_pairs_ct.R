library(stringr)
setwd("~/DataFantom")

#Expression matrices
load("Data/fantom.prom.tpm.RData")
load("Data/fantom.enh.tpm.RData")
load("Data/fantom.sample.annot.RData")

cell_types = colnames(Me_fpkm)
ct = commandArgs(trailingOnly=TRUE)[1]

mode = commandArgs(trailingOnly=TRUE)[2]

ct = "CNhs13466"
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

#cts <- sample(rownames(unique_sample_annot), size = 5)

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
build_dataset <- function (ct, folder){
  
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
    pos=  pos[pos$Enhancer %in% controled_enhs, c("Enhancer", "Gene", "Distance", "Interaction")]
    neg$interaction = 0
    colnames(neg) <-c("Enhancer", "Gene", "Distance", "Interaction")
    neg <- neg[c("Enhancer", "Gene", "Distance", "Interaction")]
  }
  par(mfrow = c(1,2))
  hist(Mg_coding_log[as.character(pos$Gene),ct], breaks = 1000)
  hist(Mg_coding_log[as.character(neg$Gene),ct], breaks = 1000)
  par(mfrow = c(1,1))
  data_train <- rbind.data.frame(pos, neg, make.row.names = TRUE )
  write.table(data_train, file = paste0(folder, "/EG_pairs_",ct,".txt"), row.names=FALSE, quote = F, sep = ';')
  return(data_train)
}

data_train <- build_dataset(ct, "Data")

length_auc_datasets <- function(){
  load("Data/Order_cell_type_pair_counts.RData")
  counts_ct = subset(counts_ct, counts_ct$pairs > 200)
  cts <- sample(rownames(counts_ct), size = 200)
  print(cts)
  print(counts_ct[cts,]$pairs)
  #cts <- read.table("cell_types_plt.txt", h = F, sep = '\t')
  #cts <- cts[-1]
  for(ct in cts){
    print(ct)
    data <- build_dataset(ct, "EG_pairs")
  }
}
print(mode)
if(mode == "plot"){
  length_auc_datasets()
}