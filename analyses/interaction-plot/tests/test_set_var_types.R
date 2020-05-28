# R script to test gene lists

library("testthat")

#get test data
script_dir <- file.path("../scripts")

#source script being tested
source(file.path(script_dir, "process_inputs.R"))

#opts list:
#	opts$include_syn
#	opts$include_noncoding
#	opts$include_nontranscribed

#inputs
#default
opts_default <- list(
	include_syn = FALSE,
	include_noncoding = FALSE,
	include_nontranscribed = FALSE
)
#include syn
opts_syn <- list(
	include_syn = TRUE,
	include_noncoding = FALSE,
	include_nontranscribed = FALSE
)
#include noncoding
opts_nonc <- list(
	include_syn = FALSE,
	include_noncoding = TRUE,
	include_nontranscribed = FALSE
)
#include nontranscribed
opts_nont <- list(
	include_syn = FALSE,
	include_noncoding = FALSE,
	include_nontranscribed = TRUE
)
#include syn and noncoding
opts_syn_nonc <- list(
	include_syn = TRUE,
	include_noncoding = TRUE,
	include_nontranscribed = FALSE
)
#include all
opts_all <- list(
	include_syn = TRUE,
	include_noncoding = TRUE,
	include_nontranscribed = TRUE
)

#expected output
#default
expect_default <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site"
)
#include syn
expect_syn <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site",
	"Silent",
	"Start_Codon_Ins",
	"Start_Codon_SNP",
	"Stop_Codon_Del",
	"De_novo_Start_InFrame",
	"De_novo_Start_OutOfFrame"
)
#include noncoding
expect_nonc <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site",
	"RNA",
	"Intron",
	"3'UTR",
	"5'UTR",
	"Splice_Region",
	"lincRNA"
)
#include nontranscribed
expect_nont <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site",
	"3'Flank",
	"5'Flank",
	"Targeted_Region"
)
#include syn and noncoding
expect_syn_nonc <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site",
	"Silent",
	"Start_Codon_Ins",
	"Start_Codon_SNP",
	"Stop_Codon_Del",
	"De_novo_Start_InFrame",
	"De_novo_Start_OutOfFrame",
	"RNA",
	"Intron",
	"3'UTR",
	"5'UTR",
	"Splice_Region",
	"lincRNA"
)
#include all
expect_all <- c(
	"Missense_Mutation",
	"Frame_Shift_Del",
	"In_Frame_Ins",
	"Frame_Shift_Ins",
	"Splice_Site",
	"Nonsense_Mutation",
	"In_Frame_Del",
	"Nonstop_Mutation",
	"Translation_Start_Site",
	"Silent",
	"Start_Codon_Ins",
	"Start_Codon_SNP",
	"Stop_Codon_Del",
	"De_novo_Start_InFrame",
	"De_novo_Start_OutOfFrame",
	"RNA",
	"Intron",
	"3'UTR",
	"5'UTR",
	"Splice_Region",
	"lincRNA",
	"3'Flank",
	"5'Flank",
	"Targeted_Region"
)

test_that ("types", {
  expect_identical(set_var_types(opts_default), expect_default)
  expect_identical(set_var_types(opts_syn), expect_syn)
  expect_identical(set_var_types(opts_nonc), expect_nonc)
  expect_identical(set_var_types(opts_nont), expect_nont)
  expect_identical(set_var_types(opts_syn_nonc), expect_syn_nonc)
  expect_identical(set_var_types(opts_all), expect_all)
})
