#!/bin/env python
# Author: Teja Koganti

# python3  01_setup_mafdb.py 
#  -m /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/inputs/temp.maf \
#   -d /Users/kogantit/Documents/git_repos/d3b-bix-analysis-toolkit/analyses/TMBanalysis/output/sample_var_db.sqlite \
#   -t /Users/kogantit/Documents/TMB/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed \
#   -c /Users/kogantit/Documents/TMB/gencode.v34.annotation.exononly.bed 
#   -b mutect2 \
#   -i ../output/xgen_intersected_cds.bed
#
# This script does the following - 
#    1. Reads the MAF file as chunks and for each chunk -
            # Calculates VAF
            # Filters based on variant_classification columm
            # Uses oybedtools to onlyc use variants overlapping target
#   2. Establishes a connection to the database and write a table names mutect2 
#
# Note: requires pandas, sqlite and pybedtools to be installed, and expects python3
# conda install -c bioconda pybedtools
# conda install -c anaconda pandas
# conda install -c anaconda sqlite


import pandas as pd
import numpy as np
import sqlite3
import pybedtools
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--maffile', required = True,
                    help = 'path to the MAF file with all samples')
parser.add_argument('-d', '--mafdatabasename', required = True,
                    help = "MAF sqlite database name")
parser.add_argument('-t', '--targetbed', required = True,
                     help = "capture target BED")
parser.add_argument('-c', '--cdsbed', required = True,
                     help = "CDS BED")
parser.add_argument('-b', '--databasetablename', required = True,
                     help = "CDS BED")
parser.add_argument('-i', '--intersectbedname', required = True,
                     help = "Name of out  BED name intersected with target and  CDS")
args = parser.parse_args()

#Second table  that will contain variants
#  that intersect with both target and coding exons
db_intersectbed_table_name=args.databasetablename+"_intersectbed"


# Function to establish a database connection; else return error message 
def create_connection(db_file):
    """ create a database connection to a SQLite database """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print("Error connecting to database")
    return(conn) 


def calculate_intersect_BED(target_BED, CDS_bed):
    intersect_bed = pybedtools.BedTool(target_BED).intersect(CDS_bed, u=True)
    return(intersect_bed)



########################################################
###########Defining MAF columns and types###############
#
#
# Datatype for each MAF column thta is being used 
maf_types = [
    ('Hugo_Symbol', 'TEXT'),
    ('Entrez_Gene_Id', 'INTEGER'),
    ('Center', 'TEXT'),
    ('NCBI_Build', 'TEXT'),
    ('Chromosome', 'TEXT'),
    ('Start_Position', 'INTEGER'),
    ('End_Position', 'INTEGER'),
    ('Variant_Classification', 'TEXT'),
    ('Variant_Type', 'TEXT'),
    ('Reference_Allele', 'TEXT'),
    ('Tumor_Seq_Allele1', 'TEXT'),
    ('Tumor_Seq_Allele2', 'TEXT'),
    ('dbSNP_RS', 'TEXT'),
    ('Tumor_Sample_Barcode', 'TEXT'),
    ('Transcript_ID', 'TEXT'),
    ('Exon_Number', 'TEXT'),
    ('t_depth', 'INTEGER'),
    ('t_ref_count', 'INTEGER'),
    ('t_alt_count', 'INTEGER'),
    ('n_depth', 'INTEGER'),
    ('n_ref_count', 'INTEGER'),
    ('n_alt_count', 'INTEGER'),
    ('Allele', 'TEXT'),
    ('flanking_bps', 'TEXT'),
    ('VAF', 'REAL')
]

# Number of lines that read in each chunk from MAF
chunksize = 1e5
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

#Variants will be filtered based on variant_classification 
#   column if they have these column values 
# Based on "Friends of Cancer Research TMB Harmonization Project paper"
#  "https://jitc.bmj.com/content/8/1/e000147#DC1"
var_class = [
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins"
]

# Matching up needed_cols and maf_types 
needed_types = [col for col in maf_types if col[0] in needed_cols]
# Joinging tuples with space in between
maf_table_def = ", ".join([" ".join(col) for col in needed_types])
#####################################################

 

    
#####################################################
########Creating a database with two tables##########
#       One where MAF variants intersect with target region
#       Second where  MAF variants intersect with target+CDS regions
#
#
# Starting a database connection     
con = create_connection(args.mafdatabasename)
# Creating a table with name as mutect2, index as an integer and 
#    maf_table_def as the column names and datatypes

#maf_create = "CREATE TABLE {}('index' INTEGER, {})".format(
#        "mutect2", maf_table_def)

maf_create = "CREATE TABLE {}('index' INTEGER, {})".format(
        args.databasetablename, maf_table_def)
maf_intersectbed_create = "CREATE TABLE {}('index' INTEGER, {})".format(
        db_intersectbed_table_name, maf_table_def)
con.execute(maf_create)
con.execute(maf_intersectbed_create)
####################################################




#####################################################
####Filling up the tables in DB with MAF variants-###
#
# Creating chunks from MAF file 
maf_chunks = pd.read_table(args.maffile,
                               na_values=["."],
                               comment="#",
                               sep="\t",
                               chunksize=chunksize)

intersected_bed = args.intersectbedname
intersected_bed = calculate_intersect_BED(args.targetbed, args.cdsbed)
#merged_intersected_bed = intersected_bed.pybedtools.BedTool.merge()
pybedtools.BedTool.to_dataframe(intersected_bed.merge()).to_csv(
    args.intersectbedname, sep="\t", header=None, index=False)

for chunk in maf_chunks:
    chunk['VAF'] = (chunk['t_alt_count'] /
                        (chunk['t_ref_count'] + chunk['t_alt_count']))
    chunk = chunk[needed_cols]
    chunk_filtered = chunk.loc[chunk.apply(lambda x: 
                        x["Variant_Classification"]  in var_class, axis=1)]
    chunk_within_exon = pybedtools.BedTool.from_dataframe(chunk_filtered).intersect(intersected_bed,
                                                                                    u=True)
    chunk_within_target = pybedtools.BedTool.from_dataframe(chunk_filtered).intersect(args.targetbed, u=True)
    chunk_target = pybedtools.BedTool.to_dataframe(chunk_within_target, names=needed_cols)
    chunk_target.to_sql(args.databasetablename, con, if_exists='append')
    chunk_intersectbed = pybedtools.BedTool.to_dataframe(chunk_within_exon, names=needed_cols)
    chunk_intersectbed.to_sql(db_intersectbed_table_name, con, if_exists='append')

con.close()
#####################################################

    
