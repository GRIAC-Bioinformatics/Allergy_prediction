# ----------------------------------------------------------------------------------------------------
## Function:    prepData()
## Input:       None
## Output:      data_list; being a list of data.frames, one for each data layer
## Description: Loads and prepares data from different sources, such that it can later be merged
# ----------------------------------------------------------------------------------------------------

# *** TUTORIAL FOR USER *** #
# - This function loads the data required for the model. As this is specific tasks, no fully generic 
# - function could be written. However, we've tried to make the steps as easy and clear
# - aspossible. 
#
# - GOAL:
# - Load all data layers/sources in variable_filter seperately in correct format, such that
# - they can be used in model evaluation
# - For this, do the following for each layer
# - 1. Load the data
# - 2. If needed, use variable_filter$<data layer name> to extract the right variables from datasets
# - 3. Transform into data.frame with correct columnames and correct data format
# 
# - FURTHERMORE:
# - After loading all data, they are stored in list object data_list. The validateData() function then
# - checks whether all objects comply with the data requirements


prepData <- function(variable_filter, ...){
  
  ####################################################################################################
  ## Load phenotype data ##
  ####################################################################################################
  
  # Format requirements:
  # object name = phenotype
  # object type = data.frame (of which typeof is "list")
  # colnames    = "Allergy"
  # rownames    = sample_ids
  # datatype    = character

  # We use the phenotype Allergy, defined as either Asthma, Rhinitis or Eczema, with IgE sensitization 
  phenotype                   <- read.csv(...)
  # If needed, transform such that phenotype is set to "yes" where sample has allergy,"no" when not
  phenotype[phenotype==1]     <- "yes"
  phenotype[phenotype==0]     <- "no"
  # phenotype               <- factor(phenotype$Allergy , levels = c("yes", "no"))
  
  ####################################################################################################
  ## Load nasal methylation data ##
  ####################################################################################################
  
  if ("nasal_methyl" %in% names(variable_filter)){
    
    # Format requirements:
    # object name = nasal_methyl
    # object type = data.frame (of which typeof is "list")
    # colnames    = methylNasal_Mval.{CpG site id} (e.g. methylNasal_Mval.cg123456)
    # rownames    = sample_ids
    # datatype    = double
    
    nasal_methyl                <- read.csv(...)
    # Extract the CpG sites that are use in the prediction model
    nasal_methyl                <- nasal_methyl[ , variable_filter$nasal_methyl]
    # If nasal methylation data is still in beta values, transform to M-values
    nasal_methyl                <- log2(nasal_methyl/(1-nasal_methyl)) 
    # To avoid confusion, give append the following to the colnames
    colnames(nasal_methyl)      <- paste0("methylNasal_Mval.", colnames(nasal_methyl))
  
  }
  
  ####################################################################################################
  ## Load blood methylation data ##
  ####################################################################################################
  
  if ("blood_methyl" %in% names(variable_filter)){
    
    # Format requirements:
    # object name = blood_methyl
    # object type = data.frame (of which typeof is "list")
    # colnames    = methylBlood_Mval.{CpG site id} (e.g. methylBlood_Mval.cg123456)
    # rownames    = sample_ids
    # datatype    = double
    
    blood_methyl                <- read.csv(...)
    # Extract the CpG sites that are use in the prediction model
    blood_methyl                <- blood_methyl[ , variable_filter$blood_methyl]
    # If blood methylation data is still in beta values, transform to M-values
    blood_methyl                <- log2(blood_methyl/(1-blood_methyl)) 
    # To avoid confusion, give append the following to the colnames
    colnames(blood_methyl)      <- paste0("methylBlood_Mval.", colnames(blood_methyl))
  }
  
  ####################################################################################################
  ## Load gender data ##
  ####################################################################################################
  
  if ("gender" %in% names(variable_filter)){
    
    # Format requirements:
    # object name = gender
    # object type = data.frame (of which typeof is "list")
    # colnames    = gendermale
    # rownames    = sample_ids
    # datatype    = integer
    
    gender                      <- read.csv(...)
    # Gender must be a integer vector where male is set to 1, transform if needed
    gender                      <- data.frame(apply(gender, 2, function(x){ifelse(x=="male",as.integer(1),as.integer(0))}))
    colnames(gender)            <- c("gendermale")
  }
    
  ####################################################################################################
  ## Load age data ##
  ####################################################################################################
  
  if ("age" %in% names(variable_filter)){
    
    # Format requirements:
    # object name = age
    # object type = data.frame (of which typeof is "list")
    # colnames    = age
    # rownames    = sample_ids
    # datatype    = double
    
    age                         <- read.csv(...)
  }
    
  # Put all data objects into data_list
  data_list <- lapply(names(variable_filter), get)
  names(data_list) <- names(variable_filter)
  
  # Validate data output 
  validateData(data_list, variable_filter)
  
  return(data_list)
}




