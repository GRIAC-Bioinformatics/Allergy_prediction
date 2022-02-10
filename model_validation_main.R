######################################################################################################
### Name:        model_validation.R
### Description: Runs a single model with fixed variables and hyperparameters
### Date:        1 May 2020
### Authors:     Merlijn van Breugel
######################################################################################################

######################################################################################################
# Main function
######################################################################################################

# ----------------------------------------------------------------------------------------------------
## Function:    main()
## Input:       None
## output:      None
## Description: Performs full model validation pipeline: 
##               1. Sources needed functions and loads packages
##               2. Load model to be validated
##               3. Loads different data layers
##               4. Merges data into suitable data format
##               5. Creating train test partitions
##               6a.(Retrains model if requested)
##               6. Performs prediction based on model with new data
##               7. Evaluates prediction
##               8. Store output
## Remarks:     Annotated with # *** TO DO *** # where the user is required to performed certain steps               
## Dependency:  All functions
# ----------------------------------------------------------------------------------------------------
main <- function(){
  
  message('Started with main function!') 
  
  setwd('C:/Users/merlijnvanbreugel/Documents/GitHub/DNA-project/2_Prediction/Model_validation/external_validation_eva_pr_3cpgmodel')
  
  # Sources all functions
  if(!require(R.utils))       {install.packages("R.utils")};       library(R.utils)
  sourceDirectory(paste0(getwd(), "/support_functions"))
  
  # Load all packages
  loadPackages()
  
  # Load model to be validated 
  model <- readRDS('validation_model.rds')
  
  # This model contains only the variables
  # CpG M-value for cg20372759
  # CpG M-value for cg01870976
  # CpG M-value for cg24224501
  # Correction factors age and gender are included, but yield zero-coefficient in glmnet model
  
  # Prepare data
  
  # *** TO DO ***#
  # - Adapt the function prepData such that it contains the variables from model$variables
  # - Whereby the data from each data layer is stored in a seperate data.frame with the columnames 
  # - according to model$variables and rownames equal to the sample_ids
  # - All layers must be stored in a list object with correct layer names 
  # - The list variable_filter can be used to extract the right variables from large dataset
  
  data_list       <- prepData(variable_filter = model$model_vars)
  
  # Merge data
  data        <- mergeData(data_list)
  
  # Get ratio between case and control
  case_control_ratio_avg = table(data$y)["yes"]/length(data$y)
  print(paste0("Ratio case vs control = ", formatC(case_control_ratio_avg, digits = 3, format = "f")))
  
  # Make predictions with model and evaluate
  model_eval  <- modelRun(model = model,
                          retrain = FALSE,
                          # x_train = data$x_model_matrix,
                          # y_train = data$y,
                          x_test  = data$x_model_matrix,
                          y_test  = data$y)
  
  print(paste0("ROC AUC = ", formatC(as.numeric(model_eval$roc_auc), digits = 3, format = "f")))
  print(paste0("PRC AUC = ", formatC(as.numeric(model_eval$prc_auc), digits = 3, format = "f")))
  
  # Create and show plots
  plots <- makePerfPlots(prc_obj = model_eval$prc_obj, 
                         roc_obj = model_eval$roc_obj, 
                         show_plots = TRUE)
  
  # Create folder for output
  dir.create(file.path(paste0(getwd(), "/output")), showWarnings = FALSE)
  
  # Write output to folder
  saveRDS(model_eval, paste0(getwd(), "/output/model_output.rds"))
  saveRDS(plots, paste0(getwd(), "/output/plots.rds"))
  
  # Also write png plot output to directory
  png(paste0(getwd(), "/output/model_plots.png"), width = 1000, height = 600)
  grid.arrange(plots$roc_plot, plots$prc_plot, nrow = 1)
  dev.off()
  
  ### NEW CODE ###
  # Was not included during previous replication
  
  # For further analysis, it would be highly valuable to also inspect performance for  specific phenotype
  #  per observations (which of Asthma, Rhinitis, Eczema & Sensitization)
  # The actual phenotype data will not need to be shared. 
  
  phenotypes_full <- read.csv('phenotypes.csv')
  
  # The following code can make performance plots and statistics for each phenotype
  # This would give us more insight into how the model behaves for different disease-specific cases
  
  # Note that the ids in this function are matched with the X matrix.
  phenotypesPredPlots(y = data$y, 
                      y_pred = model_eval$pred_probs, 
                      phenotypes = phenotypes_full, 
                      modeling_sample_ids = rownames(data_list$x_model_matrix), 
                      output_dir = paste0(getwd(), "/output/phenotype_specific/"), 
                      plot_name_suffix=".png")
  
  
  # Small plotting function, that creates density plots and boxplot of the used CpG sites in the model
  # Having this information allows us to better understand and explain the peformance of the replicated model
  cpg_names <- c("cg24224501","cg20372759","cg01870976")
  createDataDistributionPlots(data=data, 
                             cpg_names=cpg_names, 
                             output_dir=paste0(getwd(), "/output/top3cpg_data_distribution/"), 
                             plot_name_suffix=".png")
  
}

if (getOption('run.main', default=TRUE)) {
  main()
}

