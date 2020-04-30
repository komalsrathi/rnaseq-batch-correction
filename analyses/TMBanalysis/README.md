### TMB analysis 

This analysis compute tumor mutation burden for different disease types 

###### Inputs needed 
      - MAF file - This file should include variants from all samples abvailable
      - BED file - Only variants within this BED file will be used to calculate
         final TMB scores 
      - Histology file - Using PBTA histology  file as template, this file 
         should have two columns `Kids_First_Biospecimen_ID` that contained the 
         sample name and `short_histology` that has the disease type that needs 
         to be used for generating TMB plots
      - Output file name - This will be  saved in `results` folder. Every line
         contains sample name, disease type and TMB  scores   

###### Usage 
    sh run_tmb.sh `maffile` `bedfile`  `histologies` 

######     