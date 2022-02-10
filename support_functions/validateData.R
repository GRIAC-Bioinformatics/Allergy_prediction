# ----------------------------------------------------------------------------------------------------
## Function:    validateData()
## Input:       data_list
## Output:      None
## Description: Checks whether data_list is in right format, otherwise an error is thrown
# ----------------------------------------------------------------------------------------------------
validateData <- function(data_list, variable_filter, ...){
  
  ## All name requirements ##
  if(!(all(names(data_list) %in% names(variable_filter)))){
    stop("Not the right names provided in data_list_list")
  }
      
  ## Check phenotype data_list ##
  phenotype_requirement = c(
    typeof(data_list$phenotype) == "list",
    colnames(data_list$phenotype) == "Allergy",
    typeof(data_list$phenotype$Allergy) == "character",
    all(unique(data_list$phenotype$Allergy) %in% c("yes", "no"))
  )
  if(!(all(phenotype_requirement))){
    stop("Phenotype does not comply with data_list requirements")
  }

  ## Check gender data_list ##
  gender_requirement = c(
    typeof(data_list$gender) == "list",
    colnames(data_list$gender) == "gendermale",
    typeof(data_list$gender$gendermale) == "integer",
    all(c(unique(data_list$gendermale)) %in% c(0,1))
  )
  if(!(all(gender_requirement))){
    stop("Gendermale does comply with data_list requirements")
  }
  
  ## Check age data_list ##
  age_requirement = c(
    typeof(data_list$age) == "list",
    colnames(data_list$age) == "age",
    typeof(data_list$age$age) == "double"
  )
  if(!(all(age_requirement))){
    stop("Age does not comply with data_list requirements")
  }
  
  ## Check nasal methylation data_list ##
  nasal_requirement = c(
    typeof(data_list$nasal_methyl) == "list",
    all(grepl("methylNasal_Mval.cg", colnames(data_list$nasal_methyl))),
    all(apply(data_list$nasal_methyl, 1, typeof)=="double")
  )
  if(!(all(nasal_requirement))){
    stop("Nasal_methyl does not comply with data_list requirements")
  }
  
  ## Check blood methylation data_list ##
  blood_requirement = c(
    typeof(data_list$blood_methyl) == "list",
    all(grepl("methylBlood_Mval.cg", colnames(data_list$blood_methyl))),
    all(apply(data_list$blood_methyl, 1, typeof)=="double")
  )
  if(!(all(blood_requirement))){
    stop("Blood_methyl does not comply with data_list requirements")
  }
  
  ## Check all rownames for overlap ##
  if (length(Reduce(intersect, lapply(data_list, rownames)))==0){
    stop("Sample ids have zero overlap")
  }else{
    print(paste0("Total overlapping number of samples = ", length(Reduce(intersect, lapply(data_list, rownames)))))
  }
  
}

