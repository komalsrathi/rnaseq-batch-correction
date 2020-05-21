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
      aes(x = label1, y = label2, color = cooccur_score)
    ) +
      geom_point(shape = 19, aes(color = cooccur_score)) +
      scale_x_discrete(
        position = "top",
        limits = xscale,
        breaks = unique(cooccur_df$label1)
      ) + # blank unused sections.
      scale_y_discrete(
        limits = yscale,
        breaks = unique(cooccur_df$label2)
      ) +
      scale_color_gradientn(
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
  ggsave(cooccur_plot, filename = plot_file)
  return(cooccur_plot)
}

plot_disease <- function (gene_disease_counts, disease_fig, plot_size,
  divergent_colors, na_color) {
  #Make gene by disease chart.
  disease_gene_df <- gene_disease_counts %>%
    dplyr::mutate(gene = factor(gene, levels = gene))
    display_diseases <- disease_gene_df %>%
    dplyr::group_by(disease) %>%
    dplyr::tally(wt = mutant_samples) %>%
    dplyr::arrange(desc(n)) %>%
    head(7) %>% # seven so we end up with 8 total for color reasons
    dplyr::pull(disease)
  disease_gene_df <- disease_gene_df %>%
    dplyr::mutate(disease_factor =
      forcats::fct_other(disease, keep = display_diseases) %>%
      forcats::fct_relevel(display_diseases)
    )
  xscale2 <- levels(disease_gene_df$gene) %>%
      c(rep("", plot_size - length(.)))
  disease_plot <- ggplot(
    disease_gene_df,
    aes(x = gene, y = mutant_samples, fill = disease_factor)) +
    geom_col(width = 0.7) +
    labs(
      x = "",
      y = "Samples with mutations",
      fill = "Diagnosis"
      ) +
    colorblindr::scale_fill_OkabeIto() +
    scale_x_discrete(
      limits = xscale2,
      breaks = disease_gene_df$gene
    ) +
    scale_y_continuous(expand = c(0, 0.5, 0.1, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.key.size = unit(1, "char"))
  #save disease plot
  ggsave(disease_plot, filename = disease_fig)
  return(disease_plot)
}

combine_plots <- function (cooccur_plot, disease_plot, combined_fig, cooccur_df) {
  #combine cooccur and disease gene plots.
  #NOT a generic function to combine any plots.
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
  #labels for y axis will be gene names, with extra spaces (at bottom) blank
  ylabels  <- cooccur_df$gene2 %>%
    as.character() %>%
    unique() %>%
    c(rep("", plot_size - length(.)), .)
  cooccur_plot2 <- cooccur_plot +
    scale_x_discrete(
      limits = xscale,
      breaks = c()
    ) +
    scale_y_discrete(
      limits = yscale,
      labels = ylabels
    ) +
    theme(
      plot.margin = unit(c(-3.5, 0, 0, 0), "char")# negative top margin to move plots together
    )
  #Move labels and themes for disease plot
  disease_plot2 <- disease_plot +
    theme(
      axis.text.x = element_text(
        angle = -90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.title.y = element_text(
        vjust = -10# keep the label close when combined
      )
    )

  #Combine plots with <patchwork>
  #Layout of the two plots will be one over the other (1 column),
  #with the upper plot 3/4 the height of the lower plot
  combined_plot <- disease_plot2 + cooccur_plot2 +
    plot_layout(ncol = 1, heights = c(3, 4)) +
    plot_annotation(tag_levels = "A") &
    theme(#add uniform labels
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )
  #save combined plot.
  ggsave(combined_plot, filename = combined_fig, width = 8, height = 14)
}
