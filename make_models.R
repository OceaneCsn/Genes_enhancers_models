library(glmnet)
library(pROC)
library(ggplot2)


train_lasso_model <- function(data_train, sequences = 'ep', distance = T){
  
  variables <- names(data_train)[which(names(data_train)!="interaction" & names(data_train)!="pairs" & names(data_train)!="chr"
                                 & names(data_train)!="position")]
  print(variables)
  if(sequences=='e'){
    variables <- variables[!grepl('P_', variables)]
  }
  if(sequences=='p'){
    variables <- variables[!grepl('E_', variables)]
  }
  if(!distance){
    variables <- variables[variables != "distance"]
  }
  #model with enhancer variables
  x  <- as.matrix(data_train[, variables])
  y <- as.matrix(data_train$interaction)
  glm <- cv.glmnet(x, y, alpha = 1, family="binomial", keep = T, nfolds = 10)
  return(list(glm, variables))
}

positivesrates = function(seuil, d, pred){
  pbin = ifelse(pred>seuil, 1, 0)
  err = d$interaction - pbin
  err2 = d$interaction + pbin
  
  TP = sum(err2==2)
  FP = sum(err == -1)
  c(TP/sum(d$interaction==1), FP/sum(d$interaction==0))
}

test_lasso_model <- function(data_test, glm, data_train, variables, vars = "nucleotides enhancers", ct, plot = T, add = F, a = NULL){
  #evaluate the models and compare them
  x_test <- as.matrix(data_test[, variables])
  y_test <- as.matrix(data_test$interaction)
  
  pred <- as.numeric(predict(glm, x_test, type = "response"))
  auc = roc(data_test$interaction, pred)$auc
  if(plot){
    load("Data/unique_annot_ct.RData")
    curve = apply(as.matrix(seq(from = 0, by = 0.005, to = 1)),1, positivesrates, d = data_test, pred = pred)
    roc <- data.frame(x = curve[2,], y = curve[1,])
    a <- ggplot(data = roc, aes(x = x, y = y))+ geom_line(aes(colour = paste0("Model with ",vars, " AUC : ", round(auc,3)))) + 
      geom_line(aes(x = seq(0,1,by = 1/(dim(roc)[1]-1)), y = seq(0,1,by =  1/(dim(roc)[1]-1)),colour = "Random")) + 
      xlab("FP rate") + ylab("TP rate")  + theme(legend.position = c(0.6,0.3), legend.justification = c(0,1)) + 
      scale_color_manual(name = "Legend",values = c( "purple", "black"))+ 
      ggtitle(paste("ROC", ct, ' (', unique_sample_annot[ct,]$lib_id, ')')) 
    return(a)
  }
  return(auc)
}


