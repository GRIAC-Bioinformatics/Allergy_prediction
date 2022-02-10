######################################################################################################
### Name:        createDataDistributionPlots.R
### Description: Creates boxplot, density plot and correlation of CpG sites in the allergy prediction model
### Date:        20 May 2020
### Authors:     Cancan Qi
######################################################################################################


######################################################################################################
# Import arguments
######################################################################################################

# -----------------------------------------------------------------------------------------------------
## Function:    readArgs()
## Input:       data, cpg_names, output_dir, plot_name_suffix, save_figures
## Output:      saved png files of data distribution and correlation
## Description: Function that makes boxplot and density plots of CpG sites in model
# -----------------------------------------------------------------------------------------------------

createDataDistributionPlots <- function(data, cpg_names=c("cg24224501","cg20372759","cg01870976"), 
                                output_dir, plot_name_suffix="_external_validation.png", save_figures=TRUE){
  
  ## or plot top 10 CpG sites from the rank product statistics
  # cpg<-c("cg20372759","cg01870976","cg24224501","cg22862094","cg16027132","cg20790648","cg24707200","cg15006973","cg22689016")
  
  ## please use beta value of DNA methylation for the plots
  ## "data" is generated from mergeData() function 
  methy <- data$x_model_matrix[,paste0('methylNasal_Mval.', cpg_names)]
  ## turn the M value back to beta value
  beta.matrix <- 2^methy/(1+2^methy)
  colnames(beta.matrix) <- cpg_names
  colnames(methy) <- cpg_names
  
  ## add allergy
  Allergy <- data$y
  df <- cbind(as.data.frame(beta.matrix),Allergy)
  df$sample<-rownames(df)
  df.plot<-gather(df, all_of(cpg_names), key = "cpg",value = "beta_value")
  
  df.plot$cpg<-as.factor(df.plot$cpg)
  
  # Also make correlation plots
  # Use simple function to do multiple times for different specifications
  createCorPlot <- function(data, name = "Correlation_heatmap_cpg_sites", dir_path){
    # Create correlation matrix for plot - based on beta and M-values values
    cormat <- round(cor(data),3) 
    melted_cormat <- melt(cormat)
    
    size <- dim(cormat)[1]
    graphical_output_file <- paste0(dir_path, "/", name, plot_name_suffix)
    png(graphical_output_file, width = 300 + size*25, height = 300 + size*25)
    ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      theme(axis.text.y = element_text(size = 12))+
      coord_fixed()
    p <- ggheatmap + 
      geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
    print(p)
    dev.off() 
  }
  
  
  if (save_figures == TRUE){
    dir_path = file.path(output_dir)
    if (!dir.exists(dir_path)){
      dir.create(dir_path, showWarnings = TRUE)
    }
    
    ## ===========
    ## boxplot 
    ## ===========
    
    ## original scale
    graphical_output_file <- paste0(dir_path, "/Cohort_3CpG_boxplot_orig_scale", plot_name_suffix)
    png(graphical_output_file, width = 600, height = 600)
    p <- ggplot(df.plot, aes(x=cpg,y=beta_value,colour=Allergy))+geom_boxplot()+
      xlab("CpG site")+ylab("DNA methylation (beta value)")
    print(p)
    dev.off()
    
    ## scale as 0-1
    graphical_output_file <- paste0(dir_path, "/Cohort_3CpG_boxplot_01_scaled", plot_name_suffix)
    png(graphical_output_file, width = 600, height = 600)
    p <- ggplot(df.plot, aes(x=cpg,y=beta_value,colour=Allergy))+geom_boxplot()+
      xlab("CpG site")+ylab("DNA methylation (beta value)")+ylim(0,1)
    print(p)
    dev.off()

    ## ===============
    ## density plot 
    ## ===============
    
    ## original scale
    graphical_output_file <- paste0(dir_path, "/Cohort_3CpG_density_orig_scale", plot_name_suffix)
    png(graphical_output_file, width = 1500, height = 600)
    p <- ggplot(df.plot, aes(x=beta_value, fill=Allergy)) +
      geom_density(alpha=0.4)+facet_wrap(~cpg)
    print(p)
    dev.off()
    
    ## original scale
    graphical_output_file <- paste0(dir_path, "/Cohort_3CpG_density_01_scaled", plot_name_suffix)
    png(graphical_output_file, width = 1500, height = 600)
    p <- ggplot(df.plot, aes(x=beta_value, fill=Allergy)) +
      geom_density(alpha=0.4)+facet_wrap(~cpg)+xlim(0,1)
    print(p)
    dev.off()
    
    ## ===============
    ## correlation plot 
    ## ===============
  
    createCorPlot(methy, name = "Correlation_heatmap_Mval_cpg_sites", dir_path)
    createCorPlot(data=data.frame(methy, "Allergy"=ifelse(data$y=='yes',1,0)), name = "Correlation_heatmap_Mval_and_allergy_cpg_sites", dir_path)
    
    createCorPlot(beta.matrix, name = "Correlation_heatmap_betaval_cpg_sites", dir_path)
    createCorPlot(data=data.frame(beta.matrix, "Allergy"=ifelse(data$y=='yes',1,0)), name = "Correlation_heatmap_betaval_and_allergy_cpg_sites", dir_path)
    
  }
  ## similar code to plot the top 10 CpG sites, only change the input beta matrix
}
 
