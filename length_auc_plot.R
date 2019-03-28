rm(list=ls())
setwd("~/DataFantom")
source("make_datasets.R")
source("make_models.R")

ct = "CNhs11827"
vars = c("nucl")

cts <- read.table("cell_types_plt.txt", h = F, sep = '\t')
cts <- cts[-1]

compare_auc_length <- function(vars){
  aucs <- c()
  sizes <- c()
  for(ct in cts[30:100]){
    print(ct)
    data <- create_dataset(ct, vars)
    sizes = c(sizes, dim(data)[1])
    sets <- train_test_split(data)
    data_train = sets[[1]]
    data_test = sets[[2]]
    model <- train_lasso_model(data_train)
    glm <- model[[1]]
    variables <- model[[2]]
    auc = test_lasso_model(data_test=data_test, glm=glm, data_train=data_train, 
                           variables = variables, vars = "nucleotides", ct=ct,
                           plot = F)
    aucs <- c(aucs, auc)
  }
  d <- data.frame(x = sizes, y = aucs)
  ggplot(data = d, aes(x = x, y = y))+ geom_point(aes(colour = "Auc values")) + 
    geom_smooth(aes(x = x, y = y, colour = "Interpolation")) +
    xlab("size") + ylab("auc")  + theme(legend.position = c(0.8,1), legend.justification = c(0,1)) + 
    scale_color_manual(name = "Legend",values = c( "purple", "pink"))+ 
    ggtitle("AUC depending on dataset size") 
}

compare_nucl_dexter <- function(){
  dexters <- c()
  dinucls <- c()
  for(ct in cts[10:20]){
    print(ct)
    #Dexter model
    data <- create_dataset(ct, "dexter")
    sets <- train_test_split(data)
    data_train = sets[[1]]
    data_test = sets[[2]]
    model <- train_lasso_model(data_train)
    glm <- model[[1]]
    variables <- model[[2]]
    dexter = test_lasso_model(data_test=data_test, glm=glm, data_train=data_train, 
                           variables = variables, vars = "nucleotides", ct=ct,
                           plot = F)
    dexters <- c(dexters, dexter)
    
    #Dinucleotide model
    data <- create_dataset(ct, "nucl")
    sets <- train_test_split(data)
    data_train = sets[[1]]
    data_test = sets[[2]]
    model <- train_lasso_model(data_train)
    glm <- model[[1]]
    variables <- model[[2]]
    dinucl = test_lasso_model(data_test=data_test, glm=glm, data_train=data_train, 
                              variables = variables, vars = "nucleotides", ct=ct,
                              plot = F)
    dinucls <- c(dinucls, dinucl)
  }
  d <- data.frame(x = dinucls, y = dexters)
  ggplot(data = d, aes(x = x, y = y))+ geom_point(aes(colour = "Auc coordinates")) +
    geom_line(aes(x = seq(0,1,by = 1/(dim(d)[1]-1)), y = seq(0,1,by =  1/(dim(d)[1]-1)),colour = "Equivalence")) + 
    xlab("dinucleotides composition") + ylab("dexter variables")  + theme(legend.position = c(0.7,0.2), legend.justification = c(0,1)) + 
    scale_color_manual(name = "Legend",values = c( "purple", "black"))+ 
    xlim(0.5, 0.9) + ylim(0.5, 0.9) +
    ggtitle("AUC for models using either dexter or dinucleotide composition") 
}

