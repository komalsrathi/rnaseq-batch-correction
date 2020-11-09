# function to create tsne plot 
tsne.plot <- function(mat, clin, var, title){
  set.seed(100) # set seed for reproducibility
  tsneOut <- Rtsne(t(log2(mat + 1)), check_duplicates = FALSE, theta = 0)
  tsneData <- data.frame(tsneOut$Y, clin)
  p <- ggplot(tsneData, aes_string('X1', 'X2',
                                   color = var)) +
    geom_jitter(size = 4, width = 0.5, height = 0.5, alpha = 0.5) +
    theme_bw() + theme_Publication2() + ggtitle(title) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.title=element_text(size=10),
          legend.text=element_text(size=8)) + guides(size = F)
  return(p)
}
