# R script to test making filtering the cooccurrence summary for signficant
#relationships

library("testthat")

script_dir <- file.path("../scripts")

#source script being tested
source(file.path(script_dir, "cooccur_functions.R"))

#p values to test
#test for 95, 90, 94, and 99 % confidence intervals
nine_five_con <- 0.05
nine_o_con <- 0.10
nine_four_con <- 0.06
nine_nine_con <- 0.01

#set up inputs
gene1 <- c(
  "H3F3A",
  "H3F3A",
  "H3F3A",
  "H3F3A",
  "H3F3A"
)
gene2 <- c(
  "TP53",
  "PTPRD",
  "ATRX",
  "PPM1D",
  "ERG"
)
p <- c(
  0.050913691, #will pass 0.1 and 0.06
  0.019893811, #will pass 0.05, 0.1, and 0.06
  0.01, #will pass all filters
  0.060285697, #will pass 0.1
  1 #won't pass any filter
)
input_df <- data.frame(gene1, gene2, p, stringsAsFactors=FALSE)

#expected results
#expected output when filtering for p <= 0.05
gene1 <- c(
  "H3F3A",
  "H3F3A"
)
gene2 <- c(
  "PTPRD",
  "ATRX"
)
p <- c(
  0.019893811,
  0.01
)
expect_nine_five_df <- data.frame(gene1, gene2, p, stringsAsFactors=FALSE)

#expected output when filtering for p <= 0.1
gene1 <- c(
  "H3F3A",
  "H3F3A",
  "H3F3A",
  "H3F3A"
)
gene2 <- c(
  "TP53",
  "PTPRD",
  "ATRX",
  "PPM1D"
)
p <- c(
  0.050913691,
  0.019893811,
  0.01,
  0.060285697
)
expect_nine_o_df <- data.frame(gene1, gene2, p, stringsAsFactors=FALSE)

#expected output when filtering for p <= 0.06
gene1 <- c(
  "H3F3A",
  "H3F3A",
  "H3F3A"
)
gene2 <- c(
  "TP53",
  "PTPRD",
  "ATRX"
)
p <- c(
  0.050913691,
  0.019893811,
  0.01
)
expect_nine_four <- data.frame(gene1, gene2, p, stringsAsFactors=FALSE)

#expected output when filtering for p <= 0.01
gene1 <- c("H3F3A")
gene2 <- c("ATRX")
p <- c(0.01)
expect_nine_nine_df <- data.frame(gene1, gene2, p, stringsAsFactors=FALSE)

test_that ("significance_filter", {
  expect_identical(sig_filter(input_df, nine_five_con), expect_nine_five_df)
  expect_identical(sig_filter(input_df, nine_o_con), expect_nine_o_df)
  expect_identical(sig_filter(input_df, nine_four_con), expect_nine_four)
  expect_identical(sig_filter(input_df, nine_nine_con), expect_nine_nine_df)
})
