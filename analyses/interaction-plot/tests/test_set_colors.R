# R script to test gene lists

library("testthat")

#get test data
script_dir <- file.path("../scripts")
data_dir <- file.path("../test_data")
palette_file <- file.path(data_dir, "test_palette.tsv")

#source script being tested
source(file.path(script_dir, "process_inputs.R"))

#expected output
expect_colors <- list(divergent_colors = c(
	"#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
	na_color = "#f1f1f1")

test_that ("colors", {
  expect_identical(set_colors(palette_file), expect_colors)
})
