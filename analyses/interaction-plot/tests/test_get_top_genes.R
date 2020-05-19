# R script to test get_top_genes

library("testthat")
library("tibble")

#get test data
script_dir <- file.path("../scripts")

#source script being tested
source(file.path(script_dir, "calculate_score.R"))

#setup inputs

#KCNAB2
#LINC01134
#DNAJC11

gene <- c(
  "CAMTA1",
  "CAMTA1",
  "PRDM16",
  "PRDM16",
  "PRDM16",
  "CAMTA1",
  "PRDM16",
  "CNOT6",
  "LINC01134",
  "LINC01134"
)
sample <- c(
  "BS_1Q524P3A",
  "BS_1Q524P3B",
  "BS_1Q524P3C",
  "BS_1Q524P3D",
  "BS_1Q524P3E",
  "BS_1Q524P3F",
  "BS_1Q524P3G",
  "BS_1Q524P3H",
  "BS_1Q524P3I",
  "BS_1Q524P3J"
)
mutations <- as.integer(c(1, 1, 1, 1, 2, 1, 1, 1, 1, 1))
input <- tibble(gene, sample, mutations)

#expected results
one_expected <- c("PRDM16")
two_expected <- c("PRDM16", "CAMTA1")
three_expected <- c("PRDM16", "CAMTA1", "LINC01134")

test_that ("top_genes", {
  expect_equal(get_top_genes(input, 1, 1), one_expected)
  expect_equal(get_top_genes(input, 1, 4), one_expected)
  #the function keeps two rows when filtering for min mutants
  #even if only one gene has the required minimum, reflect this in the test
  expect_equal(get_top_genes(input, 2, 4), two_expected)
  expect_equal(get_top_genes(input, 2, 3), two_expected)
  expect_equal(get_top_genes(input, 3, 3), two_expected)
  expect_equal(get_top_genes(input, 3, 2), three_expected)
})
