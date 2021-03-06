# Author: Komal S. Rathi
# Function: Collapse to unique gene symbols

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "install_pkgs.R"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "collapse_rnaseq.R"))
outdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "output/")

# parameters
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "Expression Matrix (RSEM TPM, FPKM or Expected counts) (.rds)"),
  make_option(c("--gene_sym"),  type = "logical",
              help = "Is gene symbol present?"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.rds)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
mat <- opt$mat
gene_sym <- opt$gene_sym
outfile <- opt$outfile

# read data
mat <- readRDS(mat)


# collapse to unique gene symbols
if(gene_sym == TRUE){
  mat.collapsed <- mat %>%
    column_to_rownames('gene_symbol')
} else {
  mat.collapsed <- collapse.rnaseq(mat)
}

# save output
print(dim(mat.collapsed))
saveRDS(mat.collapsed, file = outfile)
