rm(list=ls())
setwd("~/DataFantom")
library(gridExtra)
source("make_datasets.R")
source("make_models.R")

ct = "CNhs10731"

data_nucl <- create_dataset(ct, "nucl", pos_type = '11', pairs_folder = "EG_pairs")
data_nucl[data_nucl$interaction==1,]$distance
data_nucl[data_nucl$interaction==0,]$distance


hist(data_nucl[data_nucl$interaction==1,]$distance, breaks = 50,
     main = paste("EP distance in positive and negative pairs", ct), col = rgb(0,1,0,0.4))

hist(data_nucl[data_nucl$interaction==0,]$distance, breaks = 50, col = rgb(1,0,0,0.4), add= T)


ee_ep_correlations <- function(data){
  
  
  data_nucl <- data
  data_nucl <- data_nucl[order(data_nucl$pairs),]
  enh_vars <-names(data_nucl)[grepl('E_', names(data_nucl))]
  prom_vars <-names(data_nucl)[grepl('P_', names(data_nucl))]
  
  
  cors_pos = c()
  cors_neg = c()
  for (i in seq(1:length(data_nucl$pairs))){
    enh <- str_split_fixed(data_nucl$pairs[i], '/', 2)[,1]
    prom <- str_split_fixed(data_nucl$pairs[i], '/', 2)[,2]
    
    if(data_nucl$interaction[i]==1){
      cors_pos <- c(cors_pos, cor(as.numeric(data_nucl[i,enh_vars]), as.numeric(data_nucl[i,prom_vars]), method = "spearman"))
    }
    else{
      cors_neg <- c(cors_neg, cor(as.numeric(data_nucl[i,enh_vars]), as.numeric(data_nucl[i,prom_vars]), method = "spearman"))
    }
  }
  hist(cors_pos, breaks = 50,
       main = paste("Nucleotide composition correlations between enhancers and promoters in", ct), col = rgb(0,1,0,0.4))
  
  hist(cors_neg, breaks = 50, col = rgb(1,0,0,0.4), add= T)
}


ee_ep_correlations(data_nucl)

sets_nucl <- train_test_split(data_nucl)
data_train_nucl = sets_nucl[[1]]
data_test_nucl = sets_nucl[[2]]
model_nucl <- train_lasso_model(data_train_nucl, sequences = 'p')
glm_nucl <- model_nucl[[1]]
variables_nucl <- model_nucl[[2]]

nucl <- test_lasso_model(data_test=data_test_nucl, glm=glm_nucl, data_train=data_train_nucl, 
                       variables = variables_nucl, vars = "nucleotides", ct=ct)


data_dexter <- create_dataset(ct, "dexter")

sets_dexter <- train_test_split_variation(data_dexter, delete_inter = F)
data_train_dexter = sets_dexter[[1]]
data_test_dexter = sets_dexter[[2]]
train_test_intersect(data_train = data_train_dexter, data_test = data_test_dexter)
model_dexter <- train_lasso_model(data_train_dexter, sequences = 'p')
glm_dexter <- model_dexter[[1]]
variables_dexter <- model_dexter[[2]]

#length(unique(str_split_fixed(data$pairs, '/', 2)[,1]))
#iLambdaMinq_ep = which(glm$lambda==glm$lambda.min)
#varselq =names(which(glm$glmnet.fit$beta[,iLambdaMinq_ep]!=0))


dexter_sorted_variation <- test_lasso_model(data_test=data_test_dexter, glm=glm_dexter, data_train=data_train_dexter, 
                 variables = variables_dexter, vars = "dexter", ct=ct)

sets_dexter <- train_test_split_variation(data_dexter)
data_train_dexter = sets_dexter[[1]]
data_test_dexter = sets_dexter[[2]]
train_test_intersect(data_train = data_train_dexter, data_test = data_test_dexter)
model_dexter <- train_lasso_model(data_train_dexter, sequences = 'p')
glm_dexter <- model_dexter[[1]]
variables_dexter <- model_dexter[[2]]

#length(unique(str_split_fixed(data$pairs, '/', 2)[,1]))
#iLambdaMinq_ep = which(glm$lambda==glm$lambda.min)
#varselq =names(which(glm$glmnet.fit$beta[,iLambdaMinq_ep]!=0))


dexter_sorted_variation_no_int <- test_lasso_model(data_test=data_test_dexter, glm=glm_dexter, data_train=data_train_dexter, 
                                            variables = variables_dexter, vars = "dexter", ct=ct)


sets_dexter <- train_test_split_order(data_dexter)
data_train_dexter = sets_dexter[[1]]
data_test_dexter = sets_dexter[[2]]
train_test_intersect(data_train = data_train_dexter, data_test = data_test_dexter)
model_dexter <- train_lasso_model(data_train_dexter, sequences = 'p')
glm_dexter <- model_dexter[[1]]
variables_dexter <- model_dexter[[2]]

#length(unique(str_split_fixed(data$pairs, '/', 2)[,1]))
#iLambdaMinq_ep = which(glm$lambda==glm$lambda.min)
#varselq =names(which(glm$glmnet.fit$beta[,iLambdaMinq_ep]!=0))


dexter_sorted <- test_lasso_model(data_test=data_test_dexter, glm=glm_dexter, data_train=data_train_dexter, 
                                            variables = variables_dexter, vars = "dexter", ct=ct)


sets_dexter <- train_test_split(data_dexter)
data_train_dexter = sets_dexter[[1]]
data_test_dexter = sets_dexter[[2]]
train_test_intersect(data_train = data_train_dexter, data_test = data_test_dexter)
model_dexter <- train_lasso_model(data_train_dexter, sequences = 'p')
glm_dexter <- model_dexter[[1]]
variables_dexter <- model_dexter[[2]]

#length(unique(str_split_fixed(data$pairs, '/', 2)[,1]))
#iLambdaMinq_ep = which(glm$lambda==glm$lambda.min)
#varselq =names(which(glm$glmnet.fit$beta[,iLambdaMinq_ep]!=0))


dexter <- test_lasso_model(data_test=data_test_dexter, glm=glm_dexter, data_train=data_train_dexter, 
                                            variables = variables_dexter, vars = "dexter", ct=ct)

data_pwm <- create_dataset(ct, "pwm")
sets_pwm <- train_test_split_order(data_pwm)
data_train_pwm = sets_pwm[[1]]
data_test_pwm = sets_pwm[[2]]
model_pwm <- train_lasso_model(data_train_pwm, sequences = 'p')
glm_pwm <- model_pwm[[1]]
variables_pwm <- model_pwm[[2]]
pwm <- test_lasso_model(data_test=data_test_pwm, glm=glm_pwm, data_train=data_train_pwm, 
                 variables = variables_pwm, vars = "pwm", ct=ct)


grid.arrange(dexter, dexter_sorted_variation_no_int, dexter_sorted_variation, nrow = 2)
                                                                                                   