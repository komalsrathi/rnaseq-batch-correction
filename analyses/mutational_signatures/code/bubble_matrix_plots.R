#' Author: C. Savonen for ALSF CCDL and Krutika Gaonkar for D3b
#'
#'
#' Given the data.frame output from `run_deconstructSigs` that has the number of 
#' mutations per Mb that belong to each sample x signature combo, make a bubble 
#' matrix plot that displays this information for all samples but summarized by 
#' histology.
#' @param sig_num_df a data.frame with number of mutations per Mb that belong to each sample x signature combo
#' @param label a character string for the title of the plot to be passed to ggplot2::ggtitle
#' @param color_palette a set of colors to use for a color palette.
#' @param color_breaks a set of numeric breaks to be used for color palette, needs to be the same length as the color palette provided. 
#' @param grouping_by in the input data.frame the column name of the grouping variable `short_histology`[Default] 
#' @param is_sample TRUE/FALSE[Default] to plot per sample(s) tile plot
#' @return A bubble matrix plot with the number of mutations per Mb and proportion of tumors with a non-zero weight for all samples summarized by histology.

bubble_matrix_plot <- function(sig_num_df, 
                               label = "none", 
                               color_palette = NA, 
                               color_breaks = NA,
                               grouping_by= "short_histology",
                               is_sample=FALSE) {
  if (is_sample==FALSE){
  # Summarize the mut_per_mb column by histology
  grouped_sig_num <- sig_num_df %>%
    dplyr::group_by(!!(as.name(grouping_by)), signature) %>%
    dplyr::summarize(
      # Calculate the proportion of tumors with a non-zero weight for each signature
      prop_tumors = sum(num_mutations > 0) / length(num_mutations),
      # Calculate the median number of mutations per Md
      med_num = median(mut_per_mb, na.rm = TRUE)
    )
  
  # Turn 0's into NAs so they are not actually plotted
  grouped_sig_num$prop_tumors[which(grouped_sig_num$prop_tumors == 0)] <-  NA
  
  # If color palette is not specified use color brewer's YlGnBu
  if (any(is.na(color_palette))) {
    color_palette <- RColorBrewer::brewer.pal(9, name = "YlGnBu")
  }
  # Make the bubble matrix plot
  ggplot2::ggplot(grouped_sig_num, ggplot2::aes(x = !!(as.name(grouping_by)), 
                                                y = forcats::fct_rev(signature), 
                                                color = med_num)) +
    ggplot2::geom_point(ggplot2::aes(size = prop_tumors)) +
    ggplot2::scale_size("Proportion of Samples", 
                        range = c(0, 4)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = -.01),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(.75)),
      legend.key.size = ggplot2::unit(.5, "cm")
    ) +
    ggplot2::scale_colour_gradientn(name = "Median Number of \n Mutations per Mb", 
                                    colors = color_palette) +
    # Make labels on top
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(label)
  }else{
    grouped_sig_num <- sig_num_df %>%
      dplyr::group_by(signature,Tumor_Sample_Barcode) %>%
      dplyr::summarize(
        # Calculate the median number of mutations per Md
        med_num = median(mut_per_mb, na.rm = TRUE)
      )
    
    # If color palette is not specified use color brewer's YlGnBu
    if (any(is.na(color_palette))) {
      color_palette <- RColorBrewer::brewer.pal(9, name = "YlGnBu")
    }
    # Make the bubble matrix plot
    ggplot2::ggplot(grouped_sig_num, ggplot2::aes(x = Tumor_Sample_Barcode, 
                                                  y = forcats::fct_rev(signature))) +
      ggplot2::geom_tile(aes(fill = med_num)) +
      ggplot2::geom_text(aes(label = round(med_num, 1))) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 55, hjust = -.01),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(.75)),
        legend.key.size = ggplot2::unit(.5, "cm")
      ) +
      ggplot2::scale_fill_gradientn(name = "Median Number of \n Mutations per Mb", 
                                      colors = color_palette) +
      # Make labels on top
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::ggtitle(label)
  }
  
  
}
