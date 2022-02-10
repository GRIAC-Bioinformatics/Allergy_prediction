# ----------------------------------------------------------------------------------------------------
## Function:    modelRun()
## Input:       model, retrain, x_train, y_train, 
##              x_test, y_test, n_top_vars
## Output:      output_list
## Description: Evaluates model for single split in train/test, firstly retrains if wanted
# ----------------------------------------------------------------------------------------------------
modelRun <- function(model, retrain = FALSE, x_train = NULL, y_train = NULL,
                     x_test = NULL, y_test = NULL, ...){
  
  if (retrain == TRUE){
    fit <- train(x = x_train, 
                 y = y_train,
                 method = model$model_name, 
                 metric = model$caret_train_metric, 
                 tuneGrid = model$grid,
                 trControl = model$trctrl,
                 preProcess = model$preProcess)
  }else{
    fit <- model$fit
  }
  
  # Make predictions
  pred_probs <-  predict(object = fit, newdata = x_test, type = 'prob')[["yes"]]
  predict_validation <- factor(ifelse(pred_probs>=0.5,"yes","no"),levels = c("yes","no"))
  
  # Get ROC and corresponding AUC
  roc <- pROC::roc # To make sure the no package masks this function
  roc_obj <- roc(y_test, pred_probs, levels = c("no","yes"))
  
  # Get PRC and corresponding AUC
  prc_obj <- list(precision.at.all.recall.levels(pred_probs, ifelse(y_test=="no",0,1)))
  prc_auc <- AUPRC(prc_obj) 
  
  output_list <- list(
    model_fit          = model$fit,
    predict_validation = predict_validation,
    pred_probs         = pred_probs,
    y_table            = table(y_test),
    roc_obj            = roc_obj,
    prc_obj            = prc_obj,
    roc_auc            = roc_obj$auc, 
    prc_auc            = prc_auc)
  
  return(output_list)
} 