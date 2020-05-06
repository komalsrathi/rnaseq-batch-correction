#!/bin/env python
# Author: Teja Koganti

# python3 03_tmbplots.py  \ 
#  -t ../output/pbta-snv-mutect2-tmbscores.target.txt \
#  -o ../output/pbta-snv-mutect2.tmb.png
#
# This script does the following - 
#  1. takes an input file that has sample name, disease type and TMB score
#  2. Using seaborn module, implement `stripplot` to generate TMB plots per disease type
#  3. Calculates the median line for each disease type
#
#Note: requires pandas, matplotlib and seaborn to be installed, and expects python3
# conda install -c anaconda pandas
# conda install -c anaconda seaborn
# conda install -c conda-forge matplotlib
# conda install -c anaconda numpy


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tmbscorefile', required = True,
                     help = "File with TMB  and short_histology columns")
parser.add_argument('-o', '--outplotname', required = True,
                     help = "File where the TMB plot should be saved")
args = parser.parse_args()


#########################################################################
#########  Preparing input file #########################################

# Reading input file into pandas 
tmb_scores_df = pd.read_csv(args.tmbscorefile, sep=",")

# Only diseas types with more than 20 samples 
count_gt_20 = tmb_scores_df.groupby("short_histology").count().iloc[:,0]>20
                                #Counting samples per disease type
count_gt_20.name = 'count'      # This will rename the header because the header
                                #is same as original input file. 

# Joining the count DF and input tmb DF
new_tmb_df = tmb_scores_df.join(count_gt_20, on='short_histology')
# Choosing only the lines that show a true under count column
tmb_scores_df = new_tmb_df[new_tmb_df['count']].iloc[:,:3]
#############################################################################


#############################################################################
######################  Creatinbg TMB plots #################################

# Getting the disease names to plot 
diseasenames = tmb_scores_df['short_histology'].unique()
# This is the length of the median line plotted 
    #under each disease plot
median_width = 0.6
# Defining figure size
plt.figure(figsize = (20,10))
#This plots all the disease types
ax = sns.stripplot(x="short_histology", y="TMB", data=tmb_scores_df.sort_values("TMB"), jitter=1, linewidth=0.5)
#This loops over every xtick along with the name and plots the median value
for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
        sample_name = text.get_text()  # "X" or "Y"
        # calculate the median value for all replicates of either X or Y
        median_val = np.median(tmb_scores_df['TMB'][tmb_scores_df['short_histology']==sample_name])
        # plot horizontal lines across the column, centered on the tick
        ax.plot([tick-median_width/2, tick+median_width/2], [median_val, median_val],
                lw=4, color='k')
 

#Taking log scale of y-axis         
plt.yscale('log')
ax.set_yticklabels(['','','.01','.1','1','10','100'])
plt.xticks(fontsize=15, rotation=90)
ax.set_title('TMB scores',fontsize= 30) # title of plot
ax.set_xlabel('Disease type',fontsize = 20) #xlabel
ax.set_ylabel('Mutations per Mb', fontsize = 20)#ylabelplt.savefig('results/pbta-snv-mutect2.TMB.png')
ax.axhline(1, ls='dotted', color="gray")
ax.axhline(2, ls='dotted', color="gray")
ax.axhline(10, ls='dotted', color="gray")
ax.axhline(100, ls='dotted', color="gray")
plt.tight_layout()
plt.savefig(args.outplotname)
#################################################################################
