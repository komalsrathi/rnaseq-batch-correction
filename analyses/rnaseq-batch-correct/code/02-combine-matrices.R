# Author: Komal S. Rathi
# Function: Combine TPM matrices from various RNAseq datasets
# TPM matrices should be collapsed to unique gene symbols

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyr))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
outdir <- outdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "output/")

# parameters
option_list <- list(
  make_option(c("--matrices"), type = "character",
              help = "Comma separated list of expression matrices to combine (RSEM TPM, FPM or Expected counts) (.RDS)"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.RDS)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
matrices <- opt$matrices
matrices <- trimws(strsplit(matrices,",")[[1]]) 
outfile <- opt$outfile
outfile <- file.path(outdir, outfile)

# function to merge all input datasets
join.all <- function(...){
  x <- lapply(..., FUN = function(x) readRDS(x))
  x <- lapply(x, FUN = function(x) x %>% rownames_to_column('gene'))# apply on each element on list and convert rownames to column
  x <- join_all(x, by = 'gene', type = 'inner')
  x <- x %>%
    column_to_rownames('gene')
  return(x)
}
combined.mat <- join.all(matrices)

# save output
saveRDS(combined.mat, file = outfile)