# R script to test gene lists

library("testthat")
library("tibble")

#get test data
script_dir <- file.path("../scripts")
data_dir <- file.path("../test_data")
maf_file <- file.path(data_dir, "test_maf.tsv")
exclude_file <- file.path(data_dir, "test_exclude_genes.txt")

#source script being tested
source(file.path(script_dir, "calculate_score.R"))

#read maf file
maf_df <- data.table::fread(maf_file, data.table = FALSE)

#setup genelist for use as input
genes <- c(
"AL590822.1",
"PRDM16",
"LINC01134",
"AL805961.2",
"AL365255.1",
"KCNAB2",
"DNAJC11",
"CAMTA1"
)

#expected results
gene <- c(
  "CAMTA1",
  "DNAJC11",
  "KCNAB2",
  "LINC01134",
  "PRDM16"
)
sample <- c(
  "BS_1Q524P3B",
  "BS_1Q524P3B",
  "BS_1Q524P3B",
  "BS_1Q524P3B",
  "BS_1Q524P3B"
)

mutations <- as.integer(c(1, 1, 1, 1, 2))

expected <- tibble(gene, sample, mutations)

test_that ("sample_count", {
  expect_equal(gene_counts(maf_df, genes), expected)
})
