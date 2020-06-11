## TMB analysis

###### Authors : Teja Koganti for D3B

This analysis computes tumor mutation burden for different disease types.
It takes a single MAF file and filters variants  based filtering strategies from [Friends of Cancer research](https://jitc.bmj.com/content/8/1/e000147#DC1). TMB is computed  based on
`(filtered variant counts* 1000000) / target BED length`

### Calculate TMB
  1. Reads the MAF file -
      - Calculates VAF
      - Filters based on variant_classification column

  2. Reads in config file and defines metadata specific fields/columns and variants  list etc.

  3. Reads in a target config file and determines `experimental strategy` from metadata file. MAF  variants  are filtered within this target file and the BED length is calculated

  4. Calculates TMB
      - Calculates TMB based on `((# of variants)*1000000) / size of BED)`
      - Prints out the samplename, TMB, counts and disease type for every sample


      `Usage`: 01_calculate_tmb_targetflexible.py [-h] -i MAF -m METADATAFILE -o
                                             OUTFILENAME -w TARGETCONFIG

          optional arguments:
          -h, --help            show this help message and exit
          -i MAF, --maf MAF     path to the MAF file
          -m METADATAFILE, --metadatafile METADATAFILE
                           path to the metadata/histology file
          -o OUTFILENAME, --outfilename OUTFILENAME
                           Out file name
          -w TARGETCONFIG, --targetconfig TARGETCONFIG
                           File with experimental strategy and path to BED file
          -c CONFIGFILE, --configfile CONFIGFILE
                           calculate_tmb.cfg.txt file with columns for disease,
                           samplename, variant types etc.  

   `Inputs`  :
    1. cohort based MAF files [here](https://s3.console.aws.amazon.com/s3/buckets/d3b-bix-dev-data-bucket/hgg-dmg-integration/openPBTA-consensus/?region=us-east-1&tab=overview)
    2. config file for target  and calculating TMB  (in inputs folder)
    3. Histology files for metadata  [here](https://github.com/d3b-center/d3b-bix-analysis-toolkit/tree/master/data/histologies)

   `Output` :
   - [output/output/pbta-snv-consensus-mutation_tmb.txt](https://github.com/d3b-center/d3b-bix-analysis-toolkit/blob/feature/tmb_code/analyses/TMBanalysis/output/pbta-snv-consensus-mutation_tmb.txt)


### Plot TMB scores

 1. Takes an input file that has sample name, count, TMB and histology type
 2. Using matplotlib module to implement cumulative distribution function plot for every disease type
 3. Uses minimum number of samples under each disease to filter out disease types  
 4. Calculates the median line for each disease type

 `Usage`: 02_cumulative_freq_TMBplot.py [-h] -t TMB_SCORES -o OUTFILENAME -s
                                   MINSAMPLESTOPLOT

    optional arguments:
      -h, --help            show this help message and exit
      -t TMB_SCORES, --tmb_scores TMB_SCORES
                      file with TMB scores
      -o OUTFILENAME, --outfilename OUTFILENAME
                      Name of the out plot, no extension
      -s MINSAMPLESTOPLOT, --minsamplestoplot MINSAMPLESTOPLOT
                      Minimum samples from each histology/disease to plot

   Cumulative  distribution function plot  here
   ![](output/pbta-snv-mutect2.CFD.TMB.png)
