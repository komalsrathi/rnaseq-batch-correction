# Mutational Signatures

### Authors: C. Savonen for ALSF CCDL and K. Gaonkar for D3b

This analysis evaluates mutational signatures of the [consensus SNV callers file](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#consensus-mutation-call).

Here the signatures from [COSMIC signatures](https://cancer.sanger.ac.uk/cosmic) or signatures from [Alexandrov et al, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23945592) are evaluated for all samples using [deconstructSigs](https://github.com/raerose01/deconstructSigs) . 

### Run script:
```
Rscript run_mutational_signature.R --help

Options:
	-i MAF, --maf=MAF
		Input maf file

	-w WGS_BED, --wgs_bed=WGS_BED
		Input WGS bed

	-x WXS_BED, --wxs_bed=WXS_BED
		Input WXS bed

	-m METADATA_DF, --metadata_df=METADATA_DF
		Input metadata

	-p PALETTES_GRADIENT, --palettes_gradient=PALETTES_GRADIENT
		Input pallete

	-d IND_SAMPLE, --ind_sample=IND_SAMPLE
		Input sample list

	-g GROUPING_BY, --grouping_by=GROUPING_BY
		Subgroup in metadata

	-s SIGNATURES, --signatures=SIGNATURES
		deconstructSigs signatures :cosmic or nature2013

	-h, --help
		Show this help message and exit

```


### Functions:
**run_deconstructSigs()**
Reads in maf file and runs deconstructSigs ( [mut.to.sigs.input()](https://github.com/raerose01/deconstructSigs#muttosigsinput) and  [whichSignatures()](https://github.com/raerose01/deconstructSigs#whichsignatures)) to get a list of weights per signature from [COSMIC signatures](https://cancer.sanger.ac.uk/cosmic) (signatures="cosmic" [Default]) for each sample.

**bubble_matrix_plot()**
Reads in data.frame from run_deconstructSigs() and plots bubble plots for a given grouping (short_histology [Default]). If is_sample=TRUE is given as parameter then a heatmap with color coded mutational signature weights are provided for each sample.

**grouped_sig_barplot()**
Reads in data.frame from run_deconstructSigs() and plots barplot for each sample in a given grouping (short_histology [Default])
