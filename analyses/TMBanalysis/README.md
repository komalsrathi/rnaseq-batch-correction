### TMB analysis 

This analysis computes tumor mutation burden for different disease types 

##### 1. Create database with variants in MAF file
- The first step reads in MAF file in chunks, intersects variants with target region and only filters variants that show `Missense_Mutation`, `Nonsense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del` and `In_Frame_Ins` in variant_classification
- Writes the  MAF variants with selected columns. This helps use groupby feature based on each sample in the next step 
   `Inputs` - `OpenPBTA-analysis/data/pbta-snv-mutect2.vep.maf`
   `Code` - [code/01_setup_mafdb.ipynb](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/code/01_setup_mafdb.ipynb)
   `Output` - [output/var_db.sqlite](https://github.com/d3b-center/d3b-bix-analysis-toolkit/tree/TMBanalysis/analyses/TMBanalysis/output)


##### 2. Calculate tumor mutation burden scores 
- This step reads in variant database, uses `GROUP BY` to calculate  number of variants under each sample. Then calculates TMB score using the formula below. 
      `TMB  = (Number of nonsynonymous  mutations/1000000)*Number of bases in target  BED`         
- This also reads in a metadata file (histology for PBTA) and matches up the disease type with the sample names 
   `Inputs` - [output/var_db.sqlite](https://github.com/d3b-center/d3b-bix-analysis-toolkit/tree/TMBanalysis/analyses/TMBanalysis/output)
   `Code` - [code/02_calculate_tmb_from_mafdb.ipynb](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/code/02_calculate_tmb_from_mafdb.ipynb)
   `Output` - [output/pbta-snv-mutect2-tmbscores.txt](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/output/pbta-snv-mutect2-tmbscores.txt)  

##### 3. Plot TMB scores 
- This step reads in the TMB scores per sample and plots a jitter plot for every  disease type
- Every disease type with at least 20 samples are plotted and a median line is drawn on each disease plot.Log scale is used for y-axis    
   `Inputs` - [output/pbta-snv-mutect2-tmbscores.txt](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/output/pbta-snv-mutect2-tmbscores.txt)
   `Code` - [code/03_tmbplots.ipynb](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/code/03_tmbplots.ipynb)
   `Output` 
   https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/output/pbta-snv-mutect2.TMB.png


