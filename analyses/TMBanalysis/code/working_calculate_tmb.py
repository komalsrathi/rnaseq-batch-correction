import subprocess
import argparse
import sys


def install_package(package, package_list):
    if not package in package_list:
        print("Installing package " + package)
        pip._internal.main(["install", package])

def get_target_dict(target_config):
        dict_to_return = {}
        with open(target_config, "r") as target_cfg:
                for line in target_cfg.readlines():
                        line=line.split()
                        #print(line)
                        dict_to_return[line[0]] = line[1]
        return(dict_to_return)


def calculate_tmb(grouped_df, meta_data, target_dict, out):
    samplename = np.unique(grouped_df["Tumor_Sample_Barcode"])[0]
    meta_data = meta_data.set_index("Kids_First_Biospecimen_ID")
    cols = list(grouped_df.columns)
    exp_strategy = meta_data.at[samplename,"experimental_strategy"]
    disease =  meta_data.at[samplename,"short_histology"]
    if exp_strategy in target_dict.keys():
        target_bed = target_dict.get(exp_strategy)
        #print(target_bed, disease, exp_strategy)
        maf_within_target = pybedtools.BedTool.from_dataframe(grouped_df).intersect(target_bed, u=True)
        mafdf_within_target = pybedtools.BedTool.to_dataframe(maf_within_target, names=cols)
        count = mafdf_within_target.shape[0]
        bed_length = calculate_bed_length(open(target_bed, "r"))
        tmb = (count* 1000000) / bed_length
        out_line = samplename+"\t"+exp_strategy+"\t"+disease+"\t"+str(count)+"\t"+str(bed_length)+"\t"+str(tmb)
        out.write(out_line+"\n")



def calculate_bed_length(in_bed):
    total_length = 0
    for line in in_bed:
        if line.split()[0] != "chrom":
            try:
                total_length += int(line.split()[2]) - int(line.split()[1])
            except IndexError:
                print("Check BED file formatting\n")
    return total_length


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
parser.add_argument("-i", "--maf", required=True, help="path to the MAF file")
parser.add_argument(
    "-m", "--metadatafile", required=True, help="path to the metadata/histology file"
)
parser.add_argument("-o", "--outfilename", required=True, help="Out file name")
parser.add_argument("-w", "--targetconfig", required=True, help="File with experimental strategy  and path to BED file")
args = parser.parse_args()


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
]
###########################################################


####Preparing MAF file #####################################
# Loading MAF file
maf_file = pd.read_table(args.maf, na_values=["."], comment="#", sep="\t")

# Filter  MAFbased on columns and also based on "variant_classification col"
maf_file = maf_file[needed_cols]
maf_file = maf_file.loc[
    maf_file.apply(lambda x: x["Variant_Classification"] in var_class, axis=1)
]

# Calculating VAF  column
maf_file["VAF"] = maf_file["t_alt_count"] / (
    maf_file["t_ref_count"] + maf_file["t_alt_count"]
)
###########################################################


########### Preparing target files ########################
target_dict = get_target_dict(args.targetconfig)
###########################################################


###  Groupby and calculate TMB ####
outfile = open(args.outfilename, "w")
outfile.write("Samplename\texperimental_strategy\tdisease\tcount\tbedlength\tTMB\n")
metadata_df = pd.read_csv(args.metadatafile, sep="\t")
grouped_maf = maf_file.groupby("Tumor_Sample_Barcode")
line = grouped_maf.apply(calculate_tmb, metadata_df, target_dict, outfile)
