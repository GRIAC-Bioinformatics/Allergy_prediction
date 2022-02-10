# ----------------------------------------------------------------------------------------------------
## Function:    makePerfPlots()
## Input:       prc_obj, roc_obj
## Output:      prc_plot, roc_plot
## Description: Creates ggplot plots of performance curves
# ----------------------------------------------------------------------------------------------------
makePerfPlots <- function(prc_obj, roc_obj, show_plots = TRUE, ...){
  
  roc_vals <- data.frame("fpr" = roc_obj$sensitivities,
                         "tpr" = roc_obj$specificities)
  
  roc_plot <- ggplot(roc_vals) +
    geom_line(aes(x = fpr, y = tpr), size = 1, colour = "blue") +
    # scale_color_manual(values = c(MIc.colours)) +
    ggtitle("ROC curve for evaluated model") + 
    scale_x_reverse(name = "False Positive Rate (1 - Specificity)",limits = c(1,0), 
                    breaks = seq(0,1,0.1), expand = c(0, 0)) + 
    scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), 
                       breaks = seq(0,1,0.1), expand = c(0, 0)) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) +
    coord_equal() + 
    annotate("text", x = 0.2, y = 0.05, vjust = 0, label = paste("AUC =",sprintf("%.3f",roc_obj$auc))) +
    theme(axis.ticks = element_line(color = "grey80"))
  
  prc_vals <- data.frame("prec" = prc_obj[[1]]$precision,
                         "rec" = prc_obj[[1]]$recall)
  
  prc_plot <- ggplot(prc_vals) +
    geom_line(aes(x = rec, y = prec), size = 1, colour = "blue") +
    # scale_color_manual(values = c(MIc.colours)) +
    ggtitle("Precision Recall curve for evaluated model") + 
    scale_x_continuous(name = "Recall",limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0, 0)) + 
    scale_y_continuous(name = "Precision", limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0, 0)) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) +
    coord_equal() + 
    annotate("text", x = 0.2, y = 0.05, vjust = 0, label = paste("AUC =",sprintf("%.3f",AUPRC(prc_obj) ))) +
    theme(axis.ticks = element_line(color = "grey80"))
  
  if (show_plots){
    grid.arrange(roc_plot, prc_plot, nrow = 1)
  }
  return(list("prc_plot" = prc_plot,
              "roc_plot" = roc_plot))
} 