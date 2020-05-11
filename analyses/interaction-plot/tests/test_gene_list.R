# R script to test gene lists

library("testthat")

#get test data
data_dir <- file.path("../test_data")
maf_file <- file.path(data_dir, "test_maf.tsv")
exclude_file <- file.path(data_dir, "test_exclude_genes")

#read maf file
maf_df <- data.table::fread(maf_file, data.table = FALSE)

#expected outputs
full_gene_list <- c(
"AL365255.1",
"AL590822.1",
"AL805961.2",
"CAMTA1",
"DNAJC11",
"KCNAB2",
"LINC01134",
"PRDM16"
)
trimmed_gene_list <- c(
"AL365255.1",
"AL590822.1",
"AL805961.2",
"CAMTA1",
"DNAJC11","LINC01134",
"KCNAB2",
"PRDM16"
)

testthat ("genelist", {
  expect_equal(get_gene_list(maf_df, exclude_file, 'LGAT'), full_gene_list),
  expect_equal(get_gene_list(maf_df, exclude_file, 'all'), trimmed_gene_list),
  expect_equal(get_gene_list(maf_df, exclude_file, 'ALL'), trimmed_gene_list),
  expect_equal(get_gene_list(maf_df, exclude_file, 'test'), full_gene_list),
})
