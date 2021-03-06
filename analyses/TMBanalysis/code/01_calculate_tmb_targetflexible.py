##  This script takes MAF file and  metadata  file with experimental strategy
#   and target config file and computes TMB  scores for every samplename
####### Author: Teja Koganti ##################
# Inputs  - 1. MAF file - Consensus should be performed prior to using as INputs
#               Must follow the MAF format. It's okay  if the first
#               commented line is absent
########### 2. metadata file - This should have  disease column and
#               experimental strategy. Sample should not be  repeated under
#               different experimental strategies
########### 3. target config file - This should have two columns,
#               one with experimental strategy in metadata  and another with the
#               corresponding BED file for processing variant counts


# python3 01_calculate_tmb_targetflexible.py \
#   -i pbta-snv-consensus.maf
#   -m pbta-histologies.tsv \
#   -w ../inputs/target_cfg \
#   -o ../output/pbta-snv-consensus-TMB_intarget.txt

import subprocess
import argparse
import sys


# This checks the packages to be installed if not already
#   installed by user
def install_package(package, package_list):
    if package not in package_list:
        print("Installing package " + package)
        pip._internal.main(["install", package])

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



# This function returns a dictionary where the keys are the experimental_strategy
#   and the values are the BED files for those target exp_strategies
def get_target_dict(target_config):
    dict_to_return = {}
    with open(target_config, "r") as target_cfg:
        for line in target_cfg.readlines():
            line = line.split()
            # print(line)
            dict_to_return[line[0]] = line[1]
    return dict_to_return


# This is the function where we calulate  the TMB, for every grouped samples,
#   the function gets  the corresponding target  BED file and calulates the
#   variants within those target files using pybedtools


def calculate_tmb(
    grouped_df, meta_data, target_dict, out, diseasecol, samplenamecol,
        targetcol
):
    samplename = np.unique(grouped_df["Tumor_Sample_Barcode"])[0]
    meta_data = meta_data.set_index(samplenamecol)
    cols = list(grouped_df.columns)
    exp_strategy = meta_data.at[samplename, targetcol]
    disease = meta_data.at[samplename, diseasecol]
    if exp_strategy in target_dict.keys():
        grouped_df["Start_Position"] = grouped_df.apply(
            lambda x: x["Start_Position"] - 1, axis=1
        )
        target_bed = target_dict.get(exp_strategy)
        maf_within_target = pybedtools.BedTool.from_dataframe(
            grouped_df).intersect(
                target_bed, u=True
        )
        mafdf_within_target = pybedtools.BedTool.to_dataframe(
            maf_within_target, names=cols
        )
        count = mafdf_within_target.shape[0]
        bed_length = calculate_bed_length(open(target_bed, "r"))
        tmb = (count * 1000000) / bed_length
        out_line = (
            samplename
            + "\t"
            + exp_strategy
            + "\t"
            + disease
            + "\t"
            + str(count)
            + "\t"
            + str(bed_length)
            + "\t"
            + str(tmb)
        )
        out.write(out_line + "\n")


# This function takes in BED file and calculates the total length of the BED
def calculate_bed_length(in_bed):
    total_length = 0
    for line in in_bed:
        if line.split()[0] != "chrom":
            try:
                total_length += int(line.split()[2]) - int(line.split()[1])
            except IndexError:
                print("Check BED file formatting\n")
    return total_length


############### REading in input files #################################

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--maf", required=True, help="path to the MAF file")
parser.add_argument("-m", "--metadatafile", required=True, help=
    "path to the metadata/histology file"
)
parser.add_argument("-o", "--outfilename", required=True, help="Out file name")
parser.add_argument(
    "-c",
    "--configfile",
    required=True,
    help="calculate_tmb.cfg.txt file with columns for disease, samplename, variant types etc.",
)
parser.add_argument(
    "-w",
    "--targetconfig",
    required=True,
    help="File with experimental strategy  and path to BED file",
)
args = parser.parse_args()



################################################################
# Defining MAF columns and metadata spcific fields from config #

needed_cols = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Tumor_Sample_Barcode",
]

with open(args.configfile) as configlines:
    for line in configlines.readlines():
        if line.startswith("Variants_to_use_for_TMB="):
            stripped_line = line.strip("Variants_to_use_for_TMB=").rstrip("\n")
            var_class = stripped_line.split(",")
        if line.startswith("disease_column="):
            disease_col = line.lstrip("disease_column").lstrip(
                "=").rstrip("\n")
        if line.startswith("samplename_column="):
            samplename_col = line.lstrip("samplename_column").lstrip(
                "=").rstrip("\n")
        if line.startswith("typeoftargetcolumn="):
            typeoftargetcol = line.lstrip("typeoftargetcolumn").lstrip(
                "=").rstrip("\n")

###########################################################


####### Preparing MAF file ##################################
# Loading MAF file
maf_file = pd.read_table(args.maf, na_values=["."], comment="#", sep="\t")
maf_file = maf_file.drop_duplicates(keep=False)


# Filter  MAFbased on columns and also based on "variant_classification col"
maf_file = maf_file[needed_cols]
maf_file = maf_file.loc[
    maf_file.apply(lambda x: x["Variant_Classification"] in var_class, axis=1)
]

###########################################################


########### Preparing target files ########################
target_dict = get_target_dict(args.targetconfig)
############################################################


############  Groupby and calculate TMB #####################
outfile = open(args.outfilename, "w")
outfile.write(
    "Samplename\texperimental_strategy\tdisease\tcount\tbedlength\tTMB\n")
metadata_df = pd.read_csv(args.metadatafile, sep="\t")
grouped_maf = maf_file.groupby("Tumor_Sample_Barcode")
line = grouped_maf.apply(
    calculate_tmb,
    metadata_df,
    target_dict,
    outfile,
    disease_col,
    samplename_col,
    typeoftargetcol,
)
###############################################################
