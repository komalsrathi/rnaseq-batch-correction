#!/bin/env python
# Author: Teja Koganti

#python3 02_calculate_tmb_from_mafdb.py \
#  -m ../inputs/pbta-histologies.tsv \
#  -c short_histology \
#  -d ../output/sample_var_db.sqlite \
#  -s Kids_First_Biospecimen_ID \
#  -b mutect2 \
#  -o ../output/sample_out_tmbscores.txt 
#  -i ../output/xgen_intersected_cds.bed 

# This script does the following - 
#  1. Establishes a connection with sqlite 
#  2. Groupby column name "Tumor_Sample_Barcode" and counts the number of var
            # for each sample from table "mutect2"
#  3. Loads metadata file(pbta-histologies for PBTA)   
#  4. Use join to match up sample names and disease types and prints the final DF
#
# Note: requires pandas and sqlite to be installed, and expects python3
# conda install -c anaconda pandas
# conda install -c anaconda sqlite


import pandas as pd
import numpy as np
import sqlite3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--metadatafile', required = True,
                    help = 'path to the metadata/histology')
parser.add_argument('-c', '--diseasecol', required = True,
                    help = "Disease column name in metadatafile")
parser.add_argument('-d', '--databasename', required = True,
                    help = "Database name with variant MAF  tables")
parser.add_argument('-s', '--samplecol', required = True,
                     help = "Sample column name in metadatafile")
#parser.add_argument('-c', '--cdsbed', required = True,
#                     help = "CDS BED")
parser.add_argument('-b', '--databasetablename', required = True,
                     help = "tablename within database with variants within target  region")
parser.add_argument('-o', '--outfilename', required = True,
                     help = "Out file name")
parser.add_argument('-i', '--intersectbedname', required = True,
                     help = "Intersected  BED between target  and CDS")
args = parser.parse_args()


db_intersectbed_table_name = args.databasetablename+"_intersectbed"
intersected_out_file = args.outfilename+"_intersected_targetandcds.bed"


target_bed_size=77462866    # This will be counted in bash script when it is fully setup 
intersected_target_bed_size = 76472464


# Establishing a connection to the MAF db 
def create_connection(db_file):
    """ create a database connection to a SQLite database """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print("Error connecting to database")
    return(conn)    



##############################################################
##############Preparing dataframes for printing output########
# Creating and empty DF with columns for sample name and TMB
sample_tmb_colnames = ["Samplenames", "TMB"]
sample_tmb_df = pd.DataFrame(columns=sample_tmb_colnames)
sample_tmb_intersectbed_df = pd.DataFrame(columns=sample_tmb_colnames)
###############################################################



################################################################
#########Filling output dataframes#############################
################################################################
#
# Starting a database connection     
con = create_connection(args.databasename)
# Using GROUPBY to count vars in each sample from table mutect2 -> for targetbed
var_db = con.execute('select Tumor_Sample_Barcode,count(Tumor_Sample_Barcode) from '+args.databasetablename+' GROUP BY  Tumor_Sample_Barcode')

# Calculating TMB based on the results from MAF DB -> for targetbed
for each_sample in var_db.fetchall():
    samplename = each_sample[0]
    #print(samplename,  (int(each_sample[1])*1000000)/target_bed_size)
    sample_tmb_df = sample_tmb_df.append({"Samplenames" : samplename ,
                                          "TMB" : (int(each_sample[1])*1000000)/target_bed_size},ignore_index=True)
    
# Using GROUPBY to count vars in each sample from table mutect2 -> for intersectbed
var_db_intersectbed = con.execute('select Tumor_Sample_Barcode,count(Tumor_Sample_Barcode) from '
                                  +db_intersectbed_table_name+' GROUP BY  Tumor_Sample_Barcode')
# Calculating TMB based on the results from MAF DB -> for targetbed
for each_sample in var_db_intersectbed.fetchall():
    samplename = each_sample[0]
    #print(samplename,  (int(each_sample[1])*1000000)/target_bed_size)
    sample_tmb_intersectbed_df = sample_tmb_intersectbed_df.append({"Samplenames" : samplename ,
                                          "TMB" : (int(each_sample[1])*1000000)/intersected_target_bed_size} ,ignore_index=True)
    
################################################################    
    
    
    
#################################################################
########## Creating final TMB outputs ###########################
# Reading in metadata file to extract disease types
#metadata = pd.read_csv(metadatafile, sep="\t", index_col=False)
metadata = pd.read_csv(args.metadatafile, sep="\t", index_col=False)


#Checking of all samples within databasetablename are in the metadata file 
if not (np.in1d(sample_tmb_df["Samplenames"].astype(str), metadata[args.samplecol].astype(str))).all():
    print("\n-----Error message: Samples not in metadata file!-----\n")
    sys.exit()
    
# Joining TMB results and metedata file and keeping only relevant columns 
final_tmb = sample_tmb_df.join(metadata.set_index(args.samplecol), 
                               on="Samplenames")[["Samplenames",  "TMB", args.diseasecol]]
final_tmb.to_csv(args.outfilename, index=False)
intersectedbed_final_tmb = sample_tmb_intersectbed_df.join(metadata.set_index(args.samplecol), 
                               on="Samplenames")[["Samplenames",  "TMB", args.diseasecol]]
intersectedbed_final_tmb.to_csv(intersected_out_file, index=False)
#################################################################

