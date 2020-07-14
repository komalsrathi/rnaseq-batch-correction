# Author: Komal S. Rathi
# Function: Create clinical file (batch information)

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
outdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "input/")

# parameters
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "Expression Matrix (RSEM TPM) (.rds)"),
  make_option(c("--clin"), type = "character",
              default = NA,
              help = "Existing clinical file (.rds)"),
  make_option(c("--id_col"), type = "character",
              default = NA,
              help = "Sample identifier column to be used (only use with --clin option)"),
  make_option(c("--lib_col"), type = "character",
              default = NA,
              help = "Library column to be used (only use with --clin option)"),
  make_option(c("--study_col"), type = "character",
              default = NA,
              help = "Study column to be used (only use with --clin option)"),
  make_option(c("--lib"), type = "character",
              default = NA,
              help = "Library Prep Method"),
  make_option(c("--study"), type = "character",
              default = NA,
              help = "Study identifier"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.rds)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
mat <- opt$mat
clin <- opt$clin
id_col <- opt$id_col
lib_col <- opt$lib_col
study_col <- opt$study_col
library_prep <- opt$lib
study <- opt$study
outfile <- opt$outfile
outfile <- file.path(outdir, outfile)

# read data
mat <- readRDS(mat)

# some qc
if(!is.na(clin) & !is.na(id_col) & !is.na(lib_col) & !is.na(study_col) & is.na(library_prep) & is.na(study)){
  print("Using existing clinical file")
  dummy <- 'clin'
  clin <- readRDS(clin)
} else if(is.na(clin) & is.na(id_col) & is.na(lib_col) & is.na(study_col) & !is.na(library_prep) & !is.na(study)){
  print("Using library and study info")
  dummy <- 'denovo'
} else {
  print("Please check input parameters")
  break
}

# create clinical file based on input parameters
if(dummy == 'denovo'){
  clin <- data.frame(sample_id = colnames(mat), study_id = study, library = library_prep)
  clin <- clin %>%
    mutate(tmp = sample_id,
           batch = paste0(study_id, '_', library)) %>%
    column_to_rownames('tmp')
}

# create clinical file using existing clinical file
if(dummy == 'clin'){
  clin <- clin %>%
    filter(!!as.name(id_col) %in% colnames(mat)) %>%
    dplyr::rename('sample_id' = sym(id_col),
                  'study_id' = sym(study_col),
                  'library' = sym(lib_col)) %>%
    mutate(batch = paste0(study_id, '_', library))  %>%
    dplyr::select(sample_id, study_id, library, batch) %>% 
    mutate(tmp = sample_id) %>%
    column_to_rownames('tmp')
}

# save output
saveRDS(clin, file = outfile)
