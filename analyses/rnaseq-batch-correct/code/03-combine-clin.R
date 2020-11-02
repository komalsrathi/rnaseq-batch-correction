# Author: Komal S. Rathi
# Function: Combine metadata from various RNAseq datasets

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyr))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
outdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "output/")

# parameters
option_list <- list(
  make_option(c("--clin"), type = "character",
              help = "Comma separated list of metadata to combine (.tsv)"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.tsv)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
clin <- opt$clin
clin <- trimws(strsplit(clin,",")[[1]]) 
outfile <- opt$outfile
outfile <- file.path(outdir, outfile)

# function to rbind all input metadata
combine.all <- function(...){
  x <- lapply(..., FUN = function(x) read.delim(x, stringsAsFactors = F, check.names = F))
  x <- do.call("rbind", x)
  return(x)
}
combined.clin <- combine.all(clin)

# save output
write.table(combined.clin, file = outfile, quote = F, sep = "\t", row.names = F)
