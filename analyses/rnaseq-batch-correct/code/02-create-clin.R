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
  make_option(c("--clin"), type = "character",
              default = NA,
              help = "Existing clinical file (.tsv)"),
  make_option(c("--cohort_col"), type = "character",
              default = NA,
              help = "cohort column to be used for subsetting clinical file"),
  make_option(c("--id_col"), type = "character",
              default = NA,
              help = "Sample identifier column to be used (only use with --clin option)"),
  make_option(c("--lib_col"), type = "character",
              default = NA,
              help = "Library column to be used (only use with --clin option)"),
  make_option(c("--study_col"), type = "character",
              default = NA,
              help = "Study column to be used (only use with --clin option)"),
  make_option(c("--mat"), type = "character",
              default = NA,
              help = "Expression Matrix (RSEM TPM) (.rds) (when --clin is not provided)"),
  make_option(c("--lib"), type = "character",
              default = NA,
              help = "Library prep method for all samples (when --clin is not provided)"),
  make_option(c("--study"), type = "character",
              default = NA,
              help = "Study identifier for all samples (when --clin is not provided)"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.tsv)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
clin <- opt$clin
cohort_col <- opt$cohort_col
id_col <- opt$id_col
lib_col <- opt$lib_col
study_col <- opt$study_col
mat <- opt$mat
library_prep <- opt$lib
study <- opt$study
outfile <- opt$outfile
outfile <- file.path(outdir, outfile)

# some qc
if(!is.na(clin) & !is.na(id_col) & !is.na(lib_col) & !is.na(study_col) & 
   is.na(mat) & is.na(library_prep) & is.na(study)){
  print("Using existing clinical file")
  dummy <- 'existing_clinical'
  clin <- read.delim(clin, stringsAsFactors = F, check.names = F)
} else if(is.na(clin) & is.na(id_col) & is.na(lib_col) & is.na(study_col) 
          & !is.na(mat) & !is.na(library_prep) & !is.na(study)){
  print("Using expression matrix, library and study info")
  dummy <- 'denovo_clinical'
} else {
  print("Please check input parameters")
  break
}

# create clinical file based on input parameters
if(dummy == 'denovo_clinical'){
  # read expression matrix
  mat <- readRDS(mat)
  clin <- data.frame(identifier = colnames(mat), study_id = study, library = library_prep)
  clin <- clin %>%
    mutate(batch = paste0(study_id, '_', library)) %>%
    dplyr::select(identifier, batch)
}

# create clinical file using existing clinical file
if(dummy == 'existing_clinical'){
  
  # subset if cohort_col is not NA
  if(!is.na(cohort_col)){
    clin <- clin %>%
      filter(!is.na(get(cohort_col)))
  }
  clin <- clin %>%
    filter(experimental_strategy == 'RNA-Seq') %>%
    mutate(batch = paste0(get(study_col), '_', get(lib_col)),
           identifier = get(id_col)) %>%
    dplyr::select(identifier, batch)
}

# save output
print(dim(clin))
write.table(clin, file = outfile, quote = F, sep = "\t", row.names = F)
