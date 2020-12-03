# Author: Komal S. Rathi
# Function: merge RSEM data

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(reshape2))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "install_pkgs.R"))

# parameters
option_list <- list(
  make_option(c("--rsem_path"), type = "character",
              help = "Path to all RSEM genes.results files"),
  make_option(c("--type"), type = "character",
              help = "TPM, FPKM or expected_count"),
  make_option(c("--outfile"), type = "character",
              help = "Output filename (.RDS)"))


# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
rsem_path <- opt$rsem_path
type <- opt$type
outfile <- opt$outfile

# function to merge rsem files
merge.res <- function(nm){
  sample_name <- gsub(".*/","",nm)
  sample_name <- gsub('[.].*','',sample_name)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# combine all rsem files
expDat <- list.files(path = rsem_path, pattern = "*.genes.results*", recursive = TRUE, full.names = T)
expMat <- lapply(expDat, FUN = function(x) merge.res(x))
expMat <- data.table::rbindlist(expMat)
expMat <- dcast(expMat, gene_id~sample_name, value.var = type)

# save output
saveRDS(expMat, file = outfile)

