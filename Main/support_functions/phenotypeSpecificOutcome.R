######################################################################################################
### Name:        phenotypeSpecificOutcome.R
### Description: Make plots of model performance for different phenotypes
### Date:        1 May 2020
### Authors:     Merlijn van Breugel
######################################################################################################

######################################################################################################
# Make prediction distributions split over different diseases
######################################################################################################

# ----------------------------------------------------------------------------------------------------
## Function:    phenotypesPredPlots()
## Input:       y, y_pred, phenotypes, modeling_sample_ids, output_dir, plot_name_suffix
## Output:      stoed plots and output_list, containing confucion matrices per disease-specific phenotype
## Description: Summarizes prediction probabilities and creates'risk groups' of individuals 
##                based on model probability outcome. Model performance per bucket is calculated
# ----------------------------------------------------------------------------------------------------
phenotypesPredPlots <- function(y, y_pred, phenotypes, modeling_sample_ids, output_dir, plot_name_suffix="_external_validation.png"){
  
  # Check whether correct column names are in phenotypes
  
  if (!all(which(c('Rhinitis', 'Asthma', 'Allergy', 'Sensitization') %in% colnames(phenotypes)) == c(1,2,3,4))){
    stop("Column names of phenotype df do not match requirement. Should contain 'Rhinitis', 'Asthma', 'Allergy', 'Sensitization'")
  }
  eval_df <- data.frame(y = factor(ifelse(y=="yes",1,0), levels = c(1,0)), 
                        probs = y_pred)
  
  # If needed, filter such that ids map correctly between phenotypes and y
  phenotypes = phenotypes[which(rownames(phenotypes) %in% modeling_sample_ids), ]
  if (dim(phenotypes)[1]==0){
    stop("Rownames/ids do not match between phenotypes and modeling sample ids from x matrix")
  }
  # Temp: add extra combinations to phenotypes
  phenotypes$noAllergy                    <- 1-phenotypes$Allergy
  phenotypes$AsthmaRhinitisAllergy        <- phenotypes$Rhinitis * phenotypes$Asthma * phenotypes$Allergy
  phenotypes$AsthmaEczemaAllergy          <- phenotypes$Eczema   * phenotypes$Asthma * phenotypes$Allergy
  phenotypes$EczemaRhinitisAllergy        <- phenotypes$Rhinitis * phenotypes$Eczema * phenotypes$Allergy
  phenotypes$AsthmaRhinitisEczemaAllergy  <- phenotypes$Rhinitis * phenotypes$Asthma * phenotypes$Eczema * phenotypes$Allergy
  phenotypes$onlyIgECcontrol              <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema == 0 & phenotypes$Sensitization == 1,  1, 0)
  
  phenotypes$anytwodiseases_plusIge       <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema >= 2 & phenotypes$Sensitization == 1,  1, 0)
  phenotypes$onlyonedisease_plusIge       <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema == 1 & phenotypes$Sensitization == 1,  1, 0)
  phenotypes$anytwodiseases_noIge         <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema >= 2 & phenotypes$Sensitization == 0,  1, 0)
  phenotypes$onlyonedisease_noIge         <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema == 1 & phenotypes$Sensitization == 0,  1, 0)

  phenotypes$onlyonedisease_noIge_pure_control <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema == 1 & phenotypes$Sensitization == 0,  1, 0)
  phenotypes$nodisease_plusIge_pure_control <- ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema == 0 & phenotypes$Sensitization == 1,  1, 0)
  
  # Set the pure control 
  purecontrol_ind <- which(ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema + phenotypes$Sensitization == 0,  1, 0)==1 | phenotypes$onlyonedisease_noIge_pure_control==1)
  phenotypes$onlyonedisease_noIge_pure_control[-purecontrol_ind] <- NA
  
  purecontrol_ind <- which(ifelse(phenotypes$Rhinitis + phenotypes$Asthma + phenotypes$Eczema + phenotypes$Sensitization == 0,  1, 0)==1 | phenotypes$nodisease_plusIge_pure_control==1)
  phenotypes$nodisease_plusIge_pure_control[-purecontrol_ind] <- NA
  
  output_list <- list()
  for (p in colnames(phenotypes)[colSums(phenotypes, na.rm=TRUE) > 0]){
    index <- which(phenotypes[,p] == 1)
    output_list[[paste0('ConfusionMatrix_', p)]] <- caret::confusionMatrix(data = factor(ifelse(eval_df$probs[index]>0.5,1,0), levels = c(1,0)), 
                                                                             reference = eval_df$y[index])$table
    
    tmp <- makePredPlots(pred_prob=eval_df$probs, 
                         phenotype=phenotypes[,p], 
                         phenotype_name=p, 
                         allergy_phenotype=phenotypes$Allergy,
                         output_dir=output_dir, 
                         save_figures=TRUE, 
                         plot_name_suffix=plot_name_suffix)
  }
  saveRDS(output_list, paste0(output_dir, "/phenotype_specific_confusion_matrices.rds"))
}

## Plotting function that can be easily adapted for multiple phenotypes
makePredPlots <- function(pred_prob, phenotype, phenotype_name, allergy_phenotype, 
                              output_dir, save_figures = TRUE, plot_name_suffix = ""){
  
  phenotype_labels <- factor(ifelse(phenotype==1, "yes", "no"), levels = c("yes", "no"))
  
  plot_df <- data.frame(pred_prob = pred_prob, 
                        phenotype = phenotype_labels,
                        allergy = factor(ifelse(allergy_phenotype==1, "yes", "no"), levels = c("yes", "no")))
  
  plot_df <- plot_df[complete.cases(plot_df),]
  # Probability density plot
  pred_prob_plot <- ggplot(plot_df,aes(x=pred_prob, fill=phenotype)) + geom_density(alpha=0.25) +
    scale_x_continuous(name = "Prediction probability") + 
    scale_y_continuous(name = "Density") +
    theme_bw() + 
    ggtitle(paste0("Prediction probabilities for replicated Elastic Net model '\n and samples with phenotype ", phenotype_name, 
                   '\n with prediction probabilities for allergy')) + 
    theme(panel.grid.minor = element_blank()) +
    guides(fill=guide_legend(title=phenotype_name)) +
    theme(axis.ticks = element_line(color = "grey80"))
  
  
  
  if (save_figures == TRUE){
    dir_path = file.path(output_dir)
    if (!dir.exists(dir_path)){
      dir.create(dir_path, showWarnings = TRUE)
    }
  
    graphical_output_file <- paste0(dir_path, "/Pred_probs_",phenotype_name, plot_name_suffix)
    png(graphical_output_file, width = 600, height = 600)
    plot(pred_prob_plot)
    dev.off()
  }
  
}