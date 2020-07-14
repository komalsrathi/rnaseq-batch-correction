# function to create density plot
density.plot <- function(mat, var, title, xlab){
  p <- ggplot(mat, aes_string(x = 'log2(value + 1)', fill = var)) + 
    geom_density(alpha = .3) +
    theme_bw() + theme_Publication2() +
    ggtitle(title) +
    xlab(xlab)
  return(p)
}
