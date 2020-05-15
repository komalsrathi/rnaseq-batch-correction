# R script to test gene lists

library("testthat")

#get test data
script_dir <- file.path("../scripts")
data_dir <- file.path("../test_data")
maf_file <- file.path(data_dir, "test_maf.tsv")
exclude_file <- file.path(data_dir, "test_exclude_genes.txt")

#source script being tested
source(file.path(script_dir, "process_inputs.R"))

#read maf file
maf_df <- data.table::fread(maf_file, data.table = FALSE)

#expected outputs
full_gene_list <- c(
"AL590822.1",
"PRDM16",
"LINC01134",
"AL805961.2",
"AL365255.1",
"KCNAB2",
"DNAJC11",
"CAMTA1"
)
trimmed_gene_list <- c(
"AL590822.1",
"PRDM16",
"AL805961.2",
"AL365255.1",
"KCNAB2",
"DNAJC11",
"CAMTA1"
)

test_that ("genelist", {
  expect_equal(get_gene_list(maf_df, exclude_file, 'LGAT'), full_gene_list)
  expect_equal(get_gene_list(maf_df, exclude_file, 'all'), trimmed_gene_list)
  expect_equal(get_gene_list(maf_df, exclude_file, 'ALL'), trimmed_gene_list)
  expect_equal(get_gene_list(maf_df, exclude_file, 'test'), full_gene_list)
})
