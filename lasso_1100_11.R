rm(list=ls())
setwd("~/DataFantom")
library(gridExtra)
source("make_datasets.R")
source("make_models.R")

#ct = "CNhs11881"
ct = commandArgs(trailingOnly=TRUE)[1]

data_1100 <- create_dataset(ct, c("nucl"), pos_type = '00', pairs_folder = "EG_pairs_00")
data_11 <- create_dataset(ct, c("nucl"), pos_type = '11', pairs_folder = "EG_pairs")

model_roc <- function(data_nucl){

  sets_nucl <- train_test_split_variation(data_nucl)
  data_train_nucl = sets_nucl[[1]]
  data_test_nucl = sets_nucl[[2]]
  train_test_intersect(data_test_nucl, data_train_nucl)
  
  model_nucl <- train_lasso_model(data_train_nucl, sequences = 'ep')
  glm_nucl <- model_nucl[[1]]
  variables_nucl <- model_nucl[[2]]
  
  nucl <- test_lasso_model(data_test=data_test_nucl, glm=glm_nucl, data_train=data_train_nucl, 
                           variables = variables_nucl, vars = "nucleotides", ct=ct)
  
  return(nucl)
}

print('Building models')

roc1100 <- model_roc(data_1100)
roc11 <- model_roc(data_11)

png(paste0("Figures/1100-11-Comparaison/rocs_", ct,'_',as.character(dim(data_1100)[1]), '_', as.character(dim(data_11)[1]), ".png"))
grid.arrange(roc1100, roc11, ncol = 2, nrow = 2)
dev.off()
