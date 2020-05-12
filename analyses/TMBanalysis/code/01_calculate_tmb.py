import subprocess
import argparse
import sys


def install_package(package, package_list):
    if not package in package_list:
        print("Installing package "+package)
        pip._internal.main(['install', package])
                         
                         

def maf_intersectbed(mafbed_df, targetbed, cols):
    mafbed_within_target = pybedtools.BedTool.from_dataframe(mafbed_df).intersect(
        targetbed,u=True)
    return(pybedtools.BedTool.to_dataframe(mafbed_within_target, names=cols))


def calculate_TMB(mafbed_df, bedsize, sample_col):
    grouped_maf = mafbed_df.groupby("Tumor_Sample_Barcode")["Chromosome"].count()
    grouped_maf = grouped_maf.reset_index()
    grouped_maf["TMB"] = (grouped_maf["Chromosome"]*1000000)/bedsize
    grouped_maf.columns = ["Tumor_Sample_Barcode", "Count", "TMB"]
    #df.rename(columns={"A": "a"
    grouped_maf = grouped_maf.rename(columns={"Tumor_Sample_Barcode" : sample_col})
    return(grouped_maf)


###################################################################
############# Checking if all packages are installed ##############

reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]

needed_packages = [
    "pandas",
    "numpy",
    "pybedtools",
]

for package in needed_packages:
    install_package(package, installed_packages)
    
##################################################################



#Importing packges 
#import modin.pandas as pd
import pandas as pd
import numpy as np
import pybedtools
import sys
import pip


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--maf_file', required = True,
                    help = 'path to the MAF file')
parser.add_argument('-m', '--metadatafile', required = True,
                    help = 'path to the metadata/histology file')
parser.add_argument('-d', '--diseasecol', required = True,
                    help = "Disease column name in metadatafile")
parser.add_argument('-s', '--samplecol', required = True,
                     help = "Sample column name in metadatafile")
parser.add_argument('-t', '--target_bed', required = True,
                     help = "target_bed")
parser.add_argument('-c', '--cds_bed', required = True,
                     help = "cds_bed")
parser.add_argument('-o', '--outfilename', required = True,
                     help = "Out file name")
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
    "In_Frame_Ins"
]

needed_cols = [
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'Tumor_Sample_Barcode',
    'Transcript_ID',
    't_depth',
    't_ref_count',
    't_alt_count',
    'n_depth',
    'n_ref_count',
    'n_alt_count',
    'Allele',
    'flanking_bps',
    'VAF']
###########################################################




########### Preparing and filtering MAF ###################
intersected_bed = pybedtools.BedTool(target_bed).intersect(cds_bed, u=True)
pybedtools.BedTool.to_dataframe(intersected_bed.merge()).to_csv(
    intersected_bed_out, sep="\t", header=None, index=False)


# Loading MAF file 
maf_file = pd.read_table(maf_file,
                               na_values=["."],
                               comment="#",
                               sep="\t")
# Calculating VAF  column
maf_file['VAF'] = (maf_file['t_alt_count'] /
                        (maf_file['t_ref_count'] + maf_file['t_alt_count']))

# Filter  MAFbased on columns and also based on "variant_classification col"
maf_file = maf_file[needed_cols]
maf_filtered = maf_file.loc[maf_file.apply(lambda x: 
                        x["Variant_Classification"]  in var_class, axis=1)]

###########################################################################




##################### Filtering MAF within target ###########################
maf_within_exon_and_target = maf_intersectbed(maf_filtered, intersected_bed, needed_cols)
tmb_scores_within_exon_and_target = calculate_TMB(maf_within_exon_and_target, 
                                                  intersected_target_bed_size, sample_col)


maf_within_target = maf_intersectbed(maf_filtered,target_bed, needed_cols)
tmb_scores_within_target = calculate_TMB(maf_within_target, 
                                         target_bed_size, sample_col)
##############################################################################




################### Mapping disease col and calculate TMB #######################

metadata = pd.read_csv(metadatafile, sep="\t", index_col=False, 
                       usecols=[sample_col, disease_col])

#Checking of all samples within db_name are in the metadata file 
if not (np.in1d(maf_within_target["Tumor_Sample_Barcode"].astype(str), 
                metadata[sample_col].astype(str))).all():
    print("\n-----Error message: Samples not in metadata file!-----\n")
    sys.exit()

final_tmb_target = tmb_scores_within_target.set_index(sample_col).join(
    metadata.set_index(sample_col))
final_tmb_target_and_cds = tmb_scores_within_exon_and_target.set_index(
    sample_col).join(metadata.set_index(sample_col))
##################################################################################











