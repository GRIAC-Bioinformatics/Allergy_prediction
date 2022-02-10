# ----------------------------------------------------------------------------------------------------
## Function:    createTrainTestData()
## Input:       x, y, n_runs, p_train_test, seed, stratified, cross_val_splits
## Output:      list(x_train_ls, y_train_ls, x_test_ls, y_test_ls, case_control_ratio_avg, split_sample)
## Description: Creates splitted datasets in x and y for training and testing. When cross_val_splits=T
##               then n_runs-cross-validation is applied. Otherwise random splits are made, which are
##               stratified if requested
# ----------------------------------------------------------------------------------------------------
createTrainTestData <- function(x, y, n_runs, p_train_test = 0.7, seed = 74830934, stratified = TRUE, 
                                cross_val_splits = TRUE, output_dir, ...){
  
  ## MvB: define seeds for splitting of train and test data
  set.seed(seed) # Needed to get same splits in each job run
  print(paste("Number of different train/test splits used:", n_runs))
  
  # Make lists to store all dat partitions in 
  x_train_ls      <- list()  
  y_train_ls      <- list()  
  x_test_ls       <- list()            
  y_test_ls       <- list()  
  
  if (cross_val_splits){
    split_sample <- createFolds(y, k = as.numeric(n_runs), returnTrain = TRUE)
  }else if(stratified){
    split_sample <- createDataPartition(y, times = n_runs, p = p_train_test)
  }else{
    split_sample <- lapply(seq(1:n_runs),function(x){
      sample(1:length(y), p_train_test*length(y))
    })
  }
  
  x_train_ls    <- lapply(split_sample,function(split){x[split,]})
  y_train_ls     <- lapply(split_sample,function(split){y[split]})
  x_test_ls        <- lapply(split_sample,function(split){x[-split,]})
  y_test_ls         <- lapply(split_sample,function(split){y[-split]})
  
  case_control_ratio_avg <- mean(unlist(lapply(y_test_ls,function(x){
    mean(ifelse(x=="no",0,1))})))
  print(paste0("Average ratio case vs control = ", formatC(case_control_ratio_avg, digits = 3, format = "f")))
  
  return(list(
    x_train_ls   = x_train_ls,  
    y_train_ls    = y_train_ls,   
    x_test_ls = x_test_ls,
    y_test_ls  = y_test_ls,
    case_control_ratio_avg = case_control_ratio_avg,
    split_sample = split_sample
  ))
}