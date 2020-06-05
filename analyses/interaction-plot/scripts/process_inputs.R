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

set_var_types <- function(opts) {
  #from the provided flags
  #return a list of variant types to include in analysis

  #set up variant types by broad category (ie synonymous / non)
  #Variant Classification with High/Moderate variant consequences from maftool
  nonsynonymous <- c(
    "Missense_Mutation",
    "Frame_Shift_Del",
    "In_Frame_Ins",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Nonsense_Mutation",
    "In_Frame_Del",
    "Nonstop_Mutation",
    "Translation_Start_Site"
  )
  #Variant Classification with Low/Modifier variant consequences
  # from maftools http://asia.ensembl.org/Help/Glossary?id=535
  synonymous <- c(
    "Silent",
    "Start_Codon_Ins",
    "Start_Codon_SNP",
    "Stop_Codon_Del",
    "De_novo_Start_InFrame",
    "De_novo_Start_OutOfFrame"
  )
  intergenic <- c("IGR")
  nontranscribed <- c(
    "3'Flank",
    "5'Flank",
    "Targeted_Region"
  )
  noncoding <- c(
    "RNA",
    "Intron",
    "3'UTR",
    "5'UTR",
    "Splice_Region",
    "lincRNA"
  )

  #determine include list based on provided options
  include <- nonsynonymous #always include nonsynonymous
  if (opts$include_syn) {
    include <- c(include, synonymous)
  }
  if (opts$include_noncoding) {
    include <- c(include, noncoding)
  }
  if (opts$include_nontranscribed) {
    include <- c(include, nontranscribed)
  }

  return(include)
}
