# R script to test filtering mutations with variant type

library("testthat")

script_dir <- file.path("../scripts")

#source script being tested
source(file.path(script_dir, "cooccur_functions.R"))

#set up input data frame
Hugo_Symbol <- c(
  "PRDM16",
  "PRDM16",
  "LINC01134",
  "KCNAB2",
  "DNAJC11"
)
Chromosome <- c(
  "chr1",
  "chr1",
  "chr1",
  "chr1",
  "chr1"
)
Start_Position <- c(
  3221772,
  3346053,
  3911292,
  6095598,
  6679956
)
Variant_Classification <- c(
  "Intron",
  "Missense_Mutation",
  "Missense_Mutation",
  "Frame_Shift_Del",
  "Silent"
)
input_df <- data.frame(Hugo_Symbol, Chromosome, Start_Position,
    Variant_Classification, stringsAsFactors=FALSE)

#variant types
#only nonsyn
include_nonsyn <- c(
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
#nonsyn and syn
include_both <- c(
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

#expected results
#only filter for nonsynonymous variants
Hugo_Symbol <- c(
  "PRDM16",
  "LINC01134",
  "KCNAB2"
)
Chromosome <- c(
  "chr1",
  "chr1",
  "chr1"
)
Start_Position <- c(
  3346053,
  3911292,
  6095598
)
Variant_Classification <- c(
  "Missense_Mutation",
  "Missense_Mutation",
  "Frame_Shift_Del"
)
expect_nonsyn <- data.frame(Hugo_Symbol, Chromosome, Start_Position,
    Variant_Classification, stringsAsFactors=FALSE)

#filter for nonsynonymous and synonymous variants
Hugo_Symbol <- c(
  "PRDM16",
  "LINC01134",
  "KCNAB2",
  "DNAJC11"
)
Chromosome <- c(
  "chr1",
  "chr1",
  "chr1",
  "chr1"
)
Start_Position <- c(
  3346053,
  3911292,
  6095598,
  6679956
)
Variant_Classification <- c(
  "Missense_Mutation",
  "Missense_Mutation",
  "Frame_Shift_Del",
  "Silent"
)
expect_both <- data.frame(Hugo_Symbol, Chromosome, Start_Position,
    Variant_Classification, stringsAsFactors=FALSE)

test_that ("mutation_filter", {
  expect_identical(filter_mutations(input_df, include_nonsyn), expect_nonsyn)
  expect_identical(filter_mutations(input_df, include_both), expect_both)
})
