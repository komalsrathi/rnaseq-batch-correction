#R script to test samples_with_disease

#get test data
script_dir <- file.path("../scripts")
data_dir <- file.path("../test_data")
meta_file <- file.path(data_dir, "test_histology.tsv")
sample_file <- file.path(data_dir, "test_sample_list.tsv")

#source script being tested
source(file.path(script_dir, "calculate_score.R"))

#read maf file
meta_df <- data.table::fread(meta_file, data.table = FALSE)
sample_df <- data.table::fread(sample_file, data.table = FALSE)

#expected results

full_sample_list <- c(
  "BS_XPZAD9PP",
  "BS_0AK4F99X",
  "BS_AHAXPFG3",
  "BS_J8EH1N7V",
  "BS_0BXY0F9N",
  "BS_D3YBV7KB",
  "BS_J9M42E4M",
  "BS_T3VJQSYM",
  "BS_2AR3AP0N",
  "BS_D5NQV7N5",
  "BS_C75Y7S77",
  "BS_23M72ABG",
  "BS_9R82A3VT",
  "BS_Y7F9E2E9"
)

hgat_sample_list <- c(
  "BS_XPZAD9PP",
  "BS_0AK4F99X",
  "BS_AHAXPFG3",
  "BS_J8EH1N7V"
)

lgat_sample_list <- c(
  "BS_D3YBV7KB",
  "BS_Y7F9E2E9"
)

#run test and compare to expected results
test_that ("sample_list", {
  expect_equal(samples_with_disease(meta_df, sample_df, 'LGAT'), lgat_sample_list)
  expect_equal(samples_with_disease(meta_df, sample_df, 'lgat'), lgat_sample_list)
  expect_equal(samples_with_disease(meta_df, sample_df, 'all'), full_sample_list)
  expect_equal(samples_with_disease(meta_df, sample_df, 'HGAT'), hgat_sample_list)
})
