# Author: Komal S. Rathi
# Function: Batch correction using sva::ComBat

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
outdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "output/")

# parameters
option_list <- list(
  make_option(c("--combined_mat"), type = "character",
              help = "Combined expression matrix with multiple batches (RSEM TPM) (.rds)"),
  make_option(c("--combined_clin"), type = "character",
              help = "Combined clinical file with multiple batches (.rds)"),
  make_option(c("--sample_id"), type = "character",
              help = "Sample identifiers matching between clinical and expression matrix"),
  make_option(c("--type"), type = "character",
              help = "Type of expression data: expected_count, FPKM or TPM"),
  make_option(c("--corrected_outfile"), type = "character",
              help = "Output filename (.RDS)"),
  make_option(c("--uncorrected_outfile"), type = "character", 
              default = NULL,
              help = "Output filename (.RDS)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
combined_mat <- opt$combined_mat
combined_clin <- opt$combined_clin
sample_id <- opt$sample_id
type <- opt$type
corrected_outfile <- opt$corrected_outfile
uncorrected_outfile <- opt$uncorrected_outfile
corrected_outfile <- file.path(outdir, corrected_outfile)

# read input data
combined_mat <- readRDS(combined_mat)
combined_clin <- read.delim(combined_clin, stringsAsFactors = F, check.names = F)
combined_clin <- combined_clin %>%
  mutate(tmp = get(sample_id)) %>%
  column_to_rownames('tmp')

# arrange sample identifiers using metadata
combined_mat <- combined_mat[,rownames(combined_clin)]

if(identical(rownames(combined_clin), colnames(combined_mat))){
  print("Matching dimensions")
} else {
  print("Check inputs")
  break
}

# uncorrected file
if(!is.null(uncorrected_outfile)){
  print("Save uncorrected file..")
  print(dim(combined_mat))
  uncorrected_outfile <- file.path(outdir, uncorrected_outfile)
  saveRDS(combined_mat, file = uncorrected_outfile)
}

# batch correct using ComBat (log2(TPM + 1))
print("Batch correct uncorrected file..")
if(type != "expected_count"){
  corrected_mat <- ComBat(dat = log2(combined_mat + 1), batch = combined_clin$batch)
  corrected_mat <- 2^(corrected_mat) # back-transform
  corrected_mat <- as.data.frame(corrected_mat)
} else {
  corrected_mat <- ComBat_seq(counts = as.matrix(combined_mat), batch = combined_clin$batch)
}

print("Save corrected file..")
print(dim(corrected_mat))
saveRDS(corrected_mat, file = corrected_outfile)
