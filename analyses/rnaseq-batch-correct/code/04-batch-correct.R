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
  make_option(c("--corrected_outfile"), type = "character",
              help = "Output filename (.RDS)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
combined_mat <- opt$combined_mat
combined_clin <- opt$combined_clin
corrected_outfile <- opt$corrected_outfile
corrected_outfile <- file.path(outdir, corrected_outfile)

# read input data
combined_mat <- readRDS(combined_mat)
combined_clin <- readRDS(combined_clin)

# batch is library + study combined 
mod = model.matrix(~as.factor(batch), data = combined_clin) 
mod0 = model.matrix(~1, data = combined_clin)

# estimate latent factors not a priority
# num.sv(combined_mat, mod, method="leek")

# arrange sample identifiers using metadata
combined_mat <- combined_mat[,rownames(combined_clin)]

if(identical(rownames(combined_clin), colnames(combined_mat))){
  print("Matching dimensions")
} else {
  print("Check inputs")
  break
}

# batch correct using ComBat (log2(TPM + 1))
corrected_mat <- ComBat(dat = log2(combined_mat + 1), batch = combined_clin$batch)
corrected_mat <- 2^(corrected_mat) # back-transform
saveRDS(corrected_mat, file = corrected_outfile)
