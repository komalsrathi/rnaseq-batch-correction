# Author: Krutika Gaonkar for D3b
# Here the signatures from COSMIC signatures is evaluated for all samples in given maf using deconstructSigs

# check if packages are available if not download 
if (!("deconstructSigs" %in% installed.packages())) {
  install.packages("deconstructSigs")
}
if (!("tidyverse" %in% installed.packages())){
  install.packages("tidyverse")
}
if (!("optparse" %in% installed.packages())){
  install.packages("optparse")
}

# load required packages
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("optparse"))


option_list <- list(
  make_option(c("-i", "--maf"),
              type = "character", default = NULL,
              help = "Input maf file"),
  make_option(c("-w", "--wgs_bed"),
              type = "character", default = NA,
              help = "Input WGS bed"),
  make_option(c("-x", "--wxs_bed"),
              type = "character", default = NA,
              help = "Input WXS bed"),
  make_option(c("-m", "--metadata_df"),
              type = "character",
              help = "Input metadata"),
  make_option(c("-p", "--palettes_gradient"),
              type = "character", default = NULL,
              help = "Input pallete"),
  make_option(c("-d", "--ind_sample"),
              type = "character", default = NULL,
              help = "Input sample list"),
  make_option(c("-g", "--grouping_by"),
              type = "character", default = NULL,
              help = "Subgroup in metadata"),
  make_option(c("-s", "--signatures"),
              type = "character", default = NULL,
              help = "deconstructSigs signatures :cosmic or nature2013")
)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

### Read in input and install packages
# rootdir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)
opt <- parse_args(OptionParser(option_list = option_list))

# read putative oncogene fusion files
maf <- read_tsv(file.path(root_dir, opt$maf)) %>% as.data.frame()

# read in wgs bed
wgs_bed <- file.path(root_dir, opt$wgs_bed)

# read in wes bed
wxs_bed <- file.path(root_dir, opt$wxs_bed)

# read in metadata
metadata_df <- read_tsv(
  file.path(root_dir, opt$metadata_df), guess_max = 10000) %>%
  dplyr::rename(Tumor_Sample_Barcode = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(
    "Tumor_Sample_Barcode",
    "experimental_strategy",
    "short_histology") %>%
  # Easier to deal with NA short histologies
  # if they are labeled something different
  dplyr::mutate(
    short_histology = as.character(tidyr::replace_na(short_histology, "none")))

# read in signature param
signatures <- opt$signatures

# read in grouping param
grouping_by <- opt$grouping_by

# source functions to make them available in this notebook; maybe be called as functions is package is built around this
# step1 run deconstructSigs from maf to signature matrix
source(file.path(root_dir, "analyses",
                 "mutational_signatures",
                 "code",
                 "run_deconstructSigs.R"))
# step2 plot matrix plots
source(file.path(root_dir, "analyses",
                 "mutational_signatures",
                 "code",
                 "bubble_matrix_plots.R"))
# step3 plot signature barplot for all samples per group
source(file.path(root_dir, "analyses",
                 "mutational_signatures",
                 "code",
                 "grouped_sigs_barplot.R"))


# Step 1. Run deconstructSigs for samples in given maf file
# run deconstructSigs and format for sample X signature dataframe
deconstructSigs_output <- run_deconstructSigs( maf = maf,
                                               wgs_bed = wgs_bed,
                                               wxs_bed = wxs_bed,
                                               metadata_df = metadata_df,
                                               grouping_by = grouping_by,
                                               signatures = signatures )
# save to file
write_tsv(deconstructSigs_output, file.path(root_dir,
                                            "analyses",
                                            "mutational_signatures",
                                            "output",
                                            paste0("deconstructSigs_",
                                                   signatures, ".tsv")))

# Set up gradient color palette for the bubble matrix plots. 
gradient_col_palette <- readr::read_tsv(
  file.path(root_dir, opt$palettes_gradient))
# Won't need NA color this time. 
gradient_col_palette <- gradient_col_palette %>%
  dplyr::filter(color_names != "na_color")


# Step 2. Bubble plot per short histology
# plot and save
label <- paste(signatures, "Signatures")
bubble_matrix_plot(deconstructSigs_output,
                   label = label,
                   color_palette = gradient_col_palette$hex_codes,
                   grouping_by = grouping_by
) + ggsave(
  file.path(root_dir, "analyses",
            "mutational_signatures",
            "output",
            paste0("bubble_matrix_", signatures, "_mutation_sig.png")),
  width = 30, height = 20, units = "cm")
# Step 3. Barplot plot per short histology
# Keep only primary tumors
ind_samples <- readr::read_tsv(file.path(root_dir, opt$ind_sample))
# filter run_deconstructSigs() output
deconstructSigs_output_primary <- deconstructSigs_output %>%
  dplyr::filter(Tumor_Sample_Barcode %in% ind_samples$Kids_First_Biospecimen_ID)
# Make grouped bar plots
lapply(unique(deconstructSigs_output$short_histology),
       grouped_sig_barplot,
       sig_num_df = deconstructSigs_output_primary,
       grouping_by = grouping_by,
       output_dir = file.path(root_dir, "analyses",
                              "mutational_signatures",
                              "output",
                              "signature_grouped_barplots"),
       label = label
)
