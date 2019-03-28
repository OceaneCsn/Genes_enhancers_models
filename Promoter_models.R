rm(list=ls())
setwd("~/DataFantom")
library(gridExtra)
source("make_datasets.R")
source("make_models.R")

ct = "CNhs11827"

data_nucl <- create_dataset(ct, "nucl")
sets_nucl <- train_test_split(data_nucl)
data_train_nucl = sets_nucl[[1]]
data_test_nucl = sets_nucl[[2]]
model_nucl <- train_lasso_model(data_train_nucl, sequences = 'p')
glm_nucl <- model_nucl[[1]]
variables_nucl <- model_nucl[[2]]

nucl <- test_lasso_model(data_test=data_test_nucl, glm=glm_nucl, data_train=data_train_nucl, 
                       variables = variables_nucl, vars = "nucleotides", ct=ct)


data_dexter <- create_dataset(ct, "dexter")
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


grid.arrange(nucl, dexter, pwm, nrow = 2)
