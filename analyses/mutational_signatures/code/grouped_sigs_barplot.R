#' Author: C. Savonen for ALSF CCDL 
#'
#'
#' Given the data.frame output from `run_deconstructSigs` where the number of mutations per Mb that belong to each sample x signature combo, plot the number of mutations per Mb for each sample as a grouped bar plot for all signatures.
#' @param sig_num_df a data.frame with number of mutations per Mb that belong to each sample x signature combo
#' @param label a character string for the title of the plot to be passed to it's png name and ggtitle
#' @param hist_groups a vector of groups to each be plotted
#' @param grouping_by in the input data.frame the column name of the grouping variable `short_histology`[Default]
#' @param output_dir where the plots should be saved.If this directory doesn' exist, will create it. Default is current directory.
#' @return A grouped barplot saved as png with the mutations per Mb for each sample from each signature


grouped_sig_barplot <- function(hist_groups, sig_num_df, output_dir = getwd(),
                                label,grouping_by="short_histology") {
  # Make the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Do this for each histology group provided in. the vector
  for (hist_group in hist_groups) {
    # Narrow the df down to just this histology's group
    histology_df <- sig_num_df %>%
      dplyr::filter(!!(as.name(grouping_by)) == hist_group, mut_per_mb > 0) 

    # Make the grouped bar plot
    histology_df %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(Tumor_Sample_Barcode, -mut_per_mb), y = mut_per_mb)) +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = signature)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1)) +
      ggplot2::ylab("Mutations per Mb") +
      ggplot2::xlab("") +
      ggplot2::ggtitle(paste(hist_group, label, "signatures"))

    # Save the plot
    ggplot2::ggsave(file.path(
      output_dir,
      paste0(
        "barplot_",
        hist_group,
        "_", label,
        "_mutation_sig.png"
      )
    ))
  }
}
