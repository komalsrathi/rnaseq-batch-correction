# Author: Komal S. Rathi
# Function: QC plots for batch correction

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "pubTheme.R"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "tsne_plots.R"))
source(file.path(root_dir, "analyses", "rnaseq-batch-correct", "util", "density_plots.R"))
plotdir <- file.path(root_dir, "analyses", "rnaseq-batch-correct", "plots/")

# parameters
option_list <- list(
  make_option(c("--uncorrected_mat"), type = "character",
              help = "Combined expression matrix with multiple batches (RSEM TPM) (.rds)"),
  make_option(c("--corrected_mat"), type = "character",
              help = "Corrected expression matrix with multiple batches (RSEM TPM) (.rds)"),
  make_option(c("--combined_clin"), type = "character",
              help = "Combined clinical file with multiple batches (.tsv)"),
  make_option(c("--sample_id"), type = "character",
              help = "Sample identifier column in clinical file matching column names in expression datasets"),
  make_option(c("--tsne_plots"), type = "character",
              help = "Summary clustering plots (.pdf)"),
  make_option(c("--density_plots"), type = "character",
              help = "Histogram of housekeeping genes (.pdf)"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
uncorrected_mat <- opt$uncorrected_mat
corrected_mat <- opt$corrected_mat
combined_clin <- opt$combined_clin
sample_id <- opt$sample_id
tsne_plots <- opt$tsne_plots
tsne_plots <- file.path(plotdir, tsne_plots)
density_plots <- opt$density_plots
density_plots <- file.path(plotdir, density_plots)

# read input data
uncorrected_mat <- readRDS(uncorrected_mat)
corrected_mat <- readRDS(corrected_mat)
combined_clin <- read.delim(combined_clin, stringsAsFactors = F, check.names = F)
combined_clin <- combined_clin %>%
  mutate(tmp = get(sample_id)) %>%
  column_to_rownames('tmp')

# arrange with clinical file
uncorrected_mat <- uncorrected_mat[,rownames(combined_clin)]
corrected_mat <- corrected_mat[,rownames(combined_clin)]

# house keeping genes
hkgenes <- c("ACTB", "TUBA1A", "TUBB", "GAPDH", "LDHA", "RPL19")

# t-SNE (uncorrected matrix)
p <- tsne.plot(mat = uncorrected_mat, clin = combined_clin, var = 'batch', title = 't-SNE before batch correction')

# house keeping genes only
uncorrected_mat.hk <- uncorrected_mat[rownames(uncorrected_mat) %in% hkgenes,]
r <- tsne.plot(mat = uncorrected_mat.hk, clin = combined_clin, var = 'batch', title = 't-SNE before batch correction\n(Housekeeping genes)')

# t-SNE (corrected matrix)
q <- tsne.plot(mat = corrected_mat, clin = combined_clin, var = 'batch', title = 't-SNE after batch correction')

# house keeping genes only
corrected_mat.hk <- corrected_mat[rownames(corrected_mat) %in% hkgenes,]
s <- tsne.plot(mat = corrected_mat.hk, clin = combined_clin, var = 'batch', title = 't-SNE after batch correction\n(Housekeeping genes)')

# save t-SNE plots
ggarrange(p, q, r, s, common.legend = T) %>%
  ggexport(filename = tsne_plots, width = 12, height = 10)

# distribution of housekeeping genes
# uncorrected  mat
uncorrected_mat.hk <- uncorrected_mat.hk %>% 
  rownames_to_column('gene') %>% 
  gather(-gene, key = !!sample_id, value = 'value') %>%
  inner_join(combined_clin, by = sample_id)
p <- density.plot(mat = uncorrected_mat.hk, 
               var = 'batch', 
               title = 'House Keeping Genes (Before ComBat correction)',
               xlab = 'log2(Uncorrected value + 1)')

# corrected mat
corrected_mat.hk <- corrected_mat.hk %>% 
  as.data.frame() %>%
  rownames_to_column('gene') %>% 
  gather(-gene, key = !!sample_id, value = 'value') %>%
  inner_join(combined_clin, by = sample_id)
q <- density.plot(mat = corrected_mat.hk, 
               var = 'batch', 
               title = 'House Keeping Genes (After ComBat correction)',
               xlab = 'log2(ComBat corrected value + 1)')

# save plots
ggarrange(p, q, ncol = 2, common.legend = T) %>%
  ggexport(filename = density_plots, width = 12, height = 6)
