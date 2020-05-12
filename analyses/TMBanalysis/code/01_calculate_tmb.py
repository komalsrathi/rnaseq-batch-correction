#!/bin/env python
# Author: Teja Koganti

# This  script does the following -
#  1. Reads the MAF file -
# 	- Calculates VAF
# 	- Filters based on variant_classification column
# 	- Uses pybedtools to only filter variants within target BED
#  2. Calculates TMB
# 	- Calculates TMB based on (# of variants)*1000000 / size of BED
# 	- Add histology type from metadata file and write to outfile

# Usage:
# python3 01_calculate_tmb.py \
# 	-i temp.maf  \
# 	-m pbta-histologies.tsv  \
# 	-d "short_histology"  \
# 	-s "Kids_First_Biospecimen_ID"  \
# 	-t xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed  \
# 	-c gencode.v34.annotation.exononly.bed  \
# 	-o OUTTEST  \
#


import subprocess
import argparse
import sys


def install_package(package, package_list):
    if not package in package_list:
        print("Installing package " + package)
        pip._internal.main(["install", package])

def maf_intersectbed(mafbed_df, targetbed, cols):
    mafbed_within_target = pybedtools.BedTool.from_dataframe(mafbed_df).intersect(
        targetbed, u=True
    )
    return pybedtools.BedTool.to_dataframe(mafbed_within_target, names=cols)

def calculate_TMB(mafbed_df, bedsize, sample_col):
    grouped_maf = mafbed_df.groupby("Tumor_Sample_Barcode")["Chromosome"].count()
    grouped_maf = grouped_maf.reset_index()
    grouped_maf["TMB"] = (grouped_maf["Chromosome"] * 1000000) / bedsize
    grouped_maf.columns = ["Tumor_Sample_Barcode", "Count", "TMB"]
    grouped_maf = grouped_maf.rename(columns={"Tumor_Sample_Barcode": sample_col})
    return grouped_maf

def calculate_bed_length(in_bed):
        total_length = 0 
        for line in in_bed:
            if(line.split()[0]!="chrom"):
                try:
                      total_length  += int(line.split()[2]) - int(line.split()[1])
                except IndexError:
                      print("Check BED file formatting\n")
        return(total_length)


###################################################################
############# Checking if all packages are installed ##############

reqs = subprocess.check_output([sys.executable, "-m", "pip", "freeze"])
installed_packages = [r.decode().split("==")[0] for r in reqs.split()]

needed_packages = [
    "pandas",
    "numpy",
    "pybedtools",
]

for package in needed_packages:
    install_package(package, installed_packages)

##################################################################


# Importing packges
# import modin.pandas as pd
import pandas as pd
import numpy as np
import pybedtools
import sys
import pip


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--maf_file", required=True, help="path to the MAF file")
parser.add_argument(
    "-m", "--metadatafile", required=True, help="path to the metadata/histology file"
)
parser.add_argument(
    "-d", "--disease_col", required=True, help="Disease column name in metadatafile"
)
parser.add_argument(
    "-s", "--sample_col", required=True, help="Sample column name in metadatafile"
)
parser.add_argument("-t", "--target_bed", required=True, help="target_bed")
parser.add_argument("-c", "--cds_bed", required=True, help="cds_bed")
parser.add_argument("-o", "--outfilename", required=True, help="Out file name")
args = parser.parse_args()


#target_bed_size = 77462866
#intersected_target_bed_size = 76472464




########################################################
###########Defining MAF columns and types###############
#
#
# Based on "Friends of Cancer Research TMB Harmonization Project paper"
# "https://jitc.bmj.com/content/8/1/e000147#DC1"

var_class = [
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
]

needed_cols = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "Tumor_Sample_Barcode",
    "Transcript_ID",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "Allele",
    "flanking_bps",
    "VAF",
]
###########################################################


###########################################################################
########### Preparing and filtering MAF ###################################
print("\nPreparing MAF dataframe... \n")
intersected_bed = pybedtools.BedTool(args.target_bed).intersect(args.cds_bed, u=True)

# Getting target  BED lengths 
target_bed_size = calculate_bed_length(open(args.target_bed, "r"))
#calculate_bed_length(open(sys.argv[1], "r"))
intersected_target_bed_size = calculate_bed_length(pybedtools.BedTool.to_dataframe(intersected_bed).set_index("chrom").to_string().split("\n"))

# Loading MAF file
maf_file = pd.read_table(args.maf_file, na_values=["."], comment="#", sep="\t")
# Calculating VAF  column
maf_file["VAF"] = maf_file["t_alt_count"] / (
    maf_file["t_ref_count"] + maf_file["t_alt_count"]
)

# Filter  MAFbased on columns and also based on "variant_classification col"
maf_file = maf_file[needed_cols]
maf_filtered = maf_file.loc[
    maf_file.apply(lambda x: x["Variant_Classification"] in var_class, axis=1)
]

###########################################################################


############################################################################
##################### Filtering MAF within target ###########################
print("\n Filtering MAF dataframe... \n")

maf_within_exon_and_target = maf_intersectbed(
    maf_filtered, intersected_bed, needed_cols
)
tmb_scores_within_exon_and_target = calculate_TMB(
    maf_within_exon_and_target, intersected_target_bed_size, args.sample_col
)

maf_within_target = maf_intersectbed(maf_filtered, args.target_bed, needed_cols)
tmb_scores_within_target = calculate_TMB(
    maf_within_target, target_bed_size, args.sample_col
)
##############################################################################


################################################################################
################### Mapping disease col and calculate TMB #######################
print("\n mapping disease names with metadata file... \n")

metadata = pd.read_csv(
    args.metadatafile,
    sep="\t",
    index_col=False,
    usecols=[args.sample_col, args.disease_col],
)

# Checking of all samples within db_name are in the metadata file
if not (
    np.in1d(
        maf_within_target["Tumor_Sample_Barcode"].astype(str),
        metadata[args.sample_col].astype(str),
    )
).all():
    print("\n-----Error message: Samples not in metadata file!-----\n")
    sys.exit()

# Matching up disease type with samplenames in TMB DF
final_tmb_target = tmb_scores_within_target.set_index(args.sample_col).join(
    metadata.set_index(args.sample_col)
)
final_tmb_target_and_cds = tmb_scores_within_exon_and_target.set_index(
    args.sample_col
).join(metadata.set_index(args.sample_col))

target_out = args.outfilename + "_withintarget.txt"
target_and_cds_out = args.outfilename + "_withintarget_and_cds.txt"

# Renaming sample column and disease type column and
# writing to output file
final_tmb_target.reset_index().rename(
    columns={args.sample_col: "Sample_name", args.disease_col: "Histology_type"}
).to_csv(target_out, index=False)
final_tmb_target_and_cds.reset_index().rename(
    columns={args.sample_col: "Sample_name", args.disease_col: "Histology_type"}
).to_csv(target_and_cds_out, index=False)

##################################################################################