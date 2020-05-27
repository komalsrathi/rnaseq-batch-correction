#' Krutika Gaonkar for D3b
#'
#'
#' Given the data.frame output from `run_deconstructSigs` that has the number of 
#' mutations per Mb that belong to each sample x signature combo, make a sample
#' matrix plot that displays this information for all samples
#' @param sig_num_df a data.frame with number of mutations per Mb that belong to each sample x signature combo
#' @param label a character string for the title of the plot to be passed to ggplot2::ggtitle
#' @param color_palette a set of colors to use for a color palette.
#' @param color_breaks a set of numeric breaks to be used for color palette, needs to be the same length as the color palette provided. 
#' @param cluster a logical param if TRUE plots cluster according to mutations per Mb associated with each signature
#' @return A per sample matrix plot with the number of mutations per Mb of tumors per signature weight for all samples or a clustering plot if cluster==TRUE

perSample_matrix_plot <- function(sig_num_df,
                               label = "none",
                               color_palette = NA,
                               color_breaks = NA,
                               cluster = FALSE) {
    sig_num_df <- sig_num_df %>%
      dplyr::group_by(signature, Tumor_Sample_Barcode) %>%
      dplyr::summarize(
        # Calculate the median number of mutations per Md
        med_num = median(mut_per_mb, na.rm = TRUE)
        )
    # If color palette is not specified use color brewer's YlGnBu
    if (any(is.na(color_palette))) {
      color_palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(9)
    }
    if ( cluster == FALSE ) {
    # Make the matrix plot
    ggplot2::ggplot(sig_num_df,
                    ggplot2::aes(x = Tumor_Sample_Barcode,
                                 y = forcats::fct_rev(signature))) +
      ggplot2::geom_tile(aes(fill = med_num)) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 55, hjust = -.01),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(.75)),
        legend.key.size = ggplot2::unit(.2, "cm"),
        legend.key.width = ggplot2::unit(.8, "cm"),
        legend.position = "top"
      ) +
      ggplot2::scale_fill_gradientn(
        name = "Median Number of \n Mutations per Mb",
        colors = color_palette) +
      # Make labels on top
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::ggtitle(label)

    } else{
     sig_num_mat <- reshape2::dcast(
      sig_num_df,
      Tumor_Sample_Barcode ~ signature,
      value.var = "med_num")
    rownames(sig_num_mat) <- sig_num_mat$Tumor_Sample_Barcode
    heatmap.2(as.matrix(sig_num_mat[,-which(colnames(sig_num_mat)=="Tumor_Sample_Barcode")]),Colv = NA,dendrogram = "row",trace = "none",col=color_palette,margins=c(12,8))
    }
}
