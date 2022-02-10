# ----------------------------------------------------------------------------------------------------
## Function:    mergeData()
## Input:       data_list; being a list of data.frames, one for each data layer
## Output:      data; list object with x, y, x_model_matrix
## Description: Merges the 
# ----------------------------------------------------------------------------------------------------
mergeData <- function(data_list, ...){
  
  df <- data_list$phenotype
  
  merge_data_names <- names(data_list)[!names(data_list) %in% "phenotype"]
  for (d in merge_data_names){
    df <- transform(merge(x = df, y = data_list[[d]], by = 'row.names', all=F), row.names=Row.names, Row.names=NULL)
  }
  df <- df[complete.cases(x), ]
  
  y <- df$Allergy
  
  x <- df[, !(colnames(df) %in% c("Allergy"))]
  
  form <- paste("~",paste(colnames(x), collapse="+"), sep="")
  x_tmp <- model.frame(x, na.action=na.pass)
  x_model_matrix<- as.matrix(model.matrix(as.formula(form), data=x_tmp))[,-1,drop=FALSE]
  
  data <- list("x"              = x,
               "x_model_matrix" = x_model_matrix,
               "y"              = y)
  return(data)
}

# phenotype <- phenotypes[, "Allergy", drop = FALSE]

