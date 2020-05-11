### TMB analysis 

This analysis computes tumor mutation burden for different disease types 

##### 1. Create database with variants in MAF file
- The first step reads in MAF file in chunks, intersects variants with target region and only filters variants that show `Missense_Mutation`, `Nonsense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del` and `In_Frame_Ins` in variant_classification
- Writes the  MAF variants with selected columns. This helps use groupby feature based on each sample in the next step 
   
   `Usage`: python3 01_setup_mafdb.py [-h] -m MAFFILE -d MAFDATABASENAME -t TARGETBED -c
                         CDSBED -b DATABASETABLENAME -i INTERSECTBEDNAME

      optional arguments:
      -h, --help            show this help message and exit
      -m MAFFILE, --maffile MAFFILE
                        path to the MAF file with all samples
      -d MAFDATABASENAME, --mafdatabasename MAFDATABASENAME
                        MAF sqlite database name
      -t TARGETBED, --targetbed TARGETBED
                        capture target BED
      -c CDSBED, --cdsbed CDSBED
                        CDS BED
      -b DATABASETABLENAME, --databasetablename DATABASETABLENAME
                        CDS BED
      -i INTERSECTBEDNAME, --intersectbedname INTERSECTBEDNAME
                        Name of out BED name intersected with target and CDS
   

   `Output` : [output/var_db.sqlite](https://github.com/d3b-center/d3b-bix-analysis-toolkit/tree/TMBanalysis/analyses/TMBanalysis/output)


##### 2. Calculate tumor mutation burden scores 
- This step reads in variant database, uses `GROUP BY` to calculate  number of variants under each sample. Then calculates TMB score using the formula below. 
      `TMB  = (Number of nonsynonymous  mutations/1000000)*Number of bases in target  BED`         
- This also reads in a metadata file (histology for PBTA) and matches up the disease type with the sample names 
   
   `Usage`: python3 02_calculate_tmb_from_mafdb.py [-h] -m METADATAFILE -c DISEASECOL -d
                                      DATABASENAME -s SAMPLECOL -b
                                      DATABASETABLENAME -o OUTFILENAME -i
                                      INTERSECTBEDNAME

      optional arguments:
      -h, --help            show this help message and exit
      -m METADATAFILE, --metadatafile METADATAFILE
                        path to the metadata/histology
      -c DISEASECOL, --diseasecol DISEASECOL
                        Disease column name in metadatafile
      -d DATABASENAME, --databasename DATABASENAME
                        Database name with variant MAF tables
      -s SAMPLECOL, --samplecol SAMPLECOL
                        Sample column name in metadatafile
      -b DATABASETABLENAME, --databasetablename DATABASETABLENAME
                        tablename within database with variants within target
                        region
      -o OUTFILENAME, --outfilename OUTFILENAME
                        Out file name
      -i INTERSECTBEDNAME, --intersectbedname INTERSECTBEDNAME
                        Intersected BED between target and CDS


   `Output` : [output/pbta-snv-mutect2-tmbscores.txt](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/TMBanalysis/analyses/TMBanalysis/output/pbta-snv-mutect2-tmbscores.target.txt)  

##### 3. Plot TMB scores 
- This step reads in the TMB scores per sample and plots a jitter plot for every  disease type
- Every disease type with at least 20 samples are plotted and a median line is drawn on each disease plot. Log scale is used for y-axis    
   `Usage` : python3 03_tmbplots.py [-h] -t TMBSCOREFILE -o OUTPLOTNAME

      optional arguments:
      -h, --help            show this help message and exit
      -t TMBSCOREFILE, --tmbscorefile TMBSCOREFILE
                        File with TMB and short_histology columns
      -o OUTPLOTNAME, --outplotname OUTPLOTNAME
                        File where the TMB plot should be saved

   `Output` :
   ![](output/pbta-snv-mutect2.TMB.png)

   Cumulative  distribution function of the same fig

   ![](output/pbta-snv-mutect2.CFD.TMB.png)


