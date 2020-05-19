#R package to generate cooccurence, gene_disease, and combined plot_size

library(ggplot2)
library(patchwork)

plot_cooccurence <- function(cooccur_df, plot_file, plot_size, divergent_colors,
  na_color){
  #take cooccur_df, generate plot, and save to file
  #create scales for consistent sizing
  #The scales will need to have plot_size elements,
  #so after getting the unique list, we concatenate on extra elements.
  #for convenience, these are just numbers 1:n
  #where n is the number of extra labels needed for the scale
  xscale <- cooccur_df$label1 %>%
    as.character() %>%
    unique() %>%
    c(1:(plot_size - length(.)))
  yscale <- cooccur_df$label2 %>%
    as.character() %>%
    unique() %>%
    #the concatenated labels need to be at the front of the Y scale,
    #since this will be at the bottom in the plot.
    c(1:(plot_size - length(.)), .)

    #make cooccur plot
    cooccur_plot <- ggplot(
      cooccur_df,
      aes(x = label1, y = label2, fill = cooccur_score)
    ) +
      geom_tile(width = 0.7, height = 0.7) +
      scale_x_discrete(
        position = "top",
        limits = xscale,
        breaks = unique(cooccur_df$label1)
      ) + # blank unused sections.
      scale_y_discrete(
        limits = yscale,
        breaks = unique(cooccur_df$label2)
      ) +
      scale_fill_gradientn(
        colors = divergent_colors,
        na.value = na_color,
        limits = c(-10, 10),
        oob = scales::squish,
      ) +
      labs(
        x = "",
        y = "",
        fill = "Co-occurence\nscore"
      ) +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        axis.text.x = element_text(
          angle = -90,
          hjust = 1,
          size = 6
        ),
        axis.text.y = element_text(size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.key.size = unit(2, "char")
      )
  print(plot_file)
  ggsave(cooccur_plot, filename = plot_file)
}
