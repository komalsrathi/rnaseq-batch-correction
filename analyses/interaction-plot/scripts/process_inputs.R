library(hash)
library(dplyr)

get_gene_list <- function(maf_df, exclude_file, disease_id) {
  #given the maf data.frame and disease id, return gene list
  #only exclude commonly mutated genes when the disease group is all
  #get initial gene list
  genes <- unique(maf_df$Hugo_Symbol)
  if (tolower(disease_id) == "all") {
    exclude_df <- data.table::fread(exclude_file, data.table = FALSE)
    genes <- genes[!genes %in% exclude_df$gene] #remove genes from exclude list
  }
  return(genes)
}

set_colors <- function(palette_file) {
  #from the pallette file, set color options
  divergent_palette <- readr::read_tsv(palette_file, col_types = readr::cols())
  divergent_colors <- divergent_palette %>%
    dplyr::filter(color_names != "na_color") %>%
    dplyr::pull(hex_codes)
  na_color <- divergent_palette %>%
    dplyr::filter(color_names == "na_color") %>%
    dplyr::pull(hex_codes)
  my_list <- list(divergent_colors = divergent_colors, na_color = na_color)
  return(my_list)
}
