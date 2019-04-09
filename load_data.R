library(stringr)
setwd("~/DataFantom")

#Expression matrices
load("fantom.prom.tpm.RData")
load("fantom.enh.tpm.RData")
load("fantom.sample.annot.RData")

cell_types = colnames(Me_fpkm)
ct = commandArgs(trailingOnly=TRUE)[1]

ct = "CNhs10729"

#fantom associations based on expression correlation
associations <- read.table("hg19_enhancer_promoter_correlations_distances_cell_type.txt", sep = '\t', h = T)


#fantom annotations to choose only TSS associated with coding genes
annot <- read.table("annotation_tss_peaks.txt", h = T, sep = '\t')
annot <- na.omit(annot)
annot_coding <- annot[which(str_split_fixed(annot$short_description, '@', 2)[,1]=="p1"),]
Mg_coding <- Mg_fpkm[rownames(Mg_fpkm) %in% annot_coding$X00Annotation,]


#choose one sample per cell type
unique_sample_annot <- sample.annot[!duplicated(sample.annot$lib_id),]
cts <- sample(rownames(unique_sample_annot), size = 5)
#unique_sample_annot[cts,]$lib_id

#builds a balanced dataset for the chosen cell type, giving E-G pairs, their distance,
#and weather their interact or not
build_dataset <- function (ct){
  
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
  # control the promoter expression so the distribution of expression is the same for the
  # promoters in the negative pairs and the positive pairs
  controled_proms = c()
  posproms = pos$Gene
  negproms = neg[neg$promoter %in% rownames(Mg_coding) & neg$promoter %in% expressed_proms,]$promoter
  negproms <- sample(negproms, size = length(negproms), replace = F)
  print(length(pos$Gene))
  print(ct)
  while(length(controled_proms) < 4/5*length(pos$Gene)){
    posprom = posproms[1]
    seuil = 0
    while(length(which(abs(Mg_coding[negproms,ct]-Mg_coding[posprom,ct]) < seuil))==0){
      seuil = seuil +1
    }
    negprom = names(which(abs(Mg_coding[negproms,ct]-Mg_coding[posprom,ct]) < seuil))[1]
    index = match(negprom, negproms)
    negproms <- negproms[-index]
    controled_proms = c(controled_proms, negprom)
    if(length(posproms)>1){
      posproms <- posproms[2:length(posproms)]
    }
    else{posproms = c()}
  }
  
  neg = neg[neg$promoter %in% controled_proms, c("enhancer", "promoter", "distance")]
  neg = neg[sample(rownames(neg), size = length(pos$Enhancer)),]
  
  hist(Mg_coding[neg$promoter,ct], breaks = 10000, xlim = c(0,1000))
  hist(Mg_coding[pos$Gene,ct], breaks = 10000, xlim = c(0,1000))
  
  neg$Interaction = 0
  colnames(neg) <-c("Enhancer", "Gene", "Distance", "Interaction")
  data_train <- rbind.data.frame(pos, neg, make.row.names = TRUE )
  return(data_train)
}

draw_hists <- function(cts){
  par(mfrow = c(2,2))
  for(ct in cts){
    data_train <- build_dataset(ct)
    hist(table(subset(data_train, data_train$Interaction == 1)$Enhancer), breaks = 20, xlim = c(0,15),
         main = paste("# genes per enhancers in", ct), col = "darksalmon")
    hist(table(subset(data_train, data_train$Interaction == 1)$Gene), breaks = 20, xlim = c(0,15), col = "cadetblue3",
         main = "# enhancers per gene, Positives")
  }
}

draw_hists_all_ct <- function(){
  par(mfrow = c(2,2))
  
  #data_train <- build_dataset(ct)
  hist(table(associations[associations$promoter %in% rownames(Mg_coding),]$enhancer), breaks = 40, xlim = c(0,30),
       main = "# genes per enhancers in", col = "darksalmon")
  hist(table(associations[associations$promoter %in% rownames(Mg_coding),]$promoter), breaks = 40, xlim = c(0,30), col = "cadetblue3",
       main = "# enhancers per gene, Positives")
  
}

#draw_hists_all_ct()
#draw_hists(cts)
data_train <- build_dataset(ct)
write.table(data_train, file = paste0("EG_pairs_",cts,".txt"), row.names=FALSE, quote = F, sep = ';')
