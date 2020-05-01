# Mutational Signatures

### Authors: C. Savonen for ALSF CCDL and K. Gaonkar for D3b

This analysis evaluates mutational signatures of the [consensus SNV callers file](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#consensus-mutation-call).

Here the signatures from [COSMIC signatures](https://cancer.sanger.ac.uk/cosmic)
and[Alexandrov et al, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23945592) are
evaluated for all samples using [deconstructSigs](https://github.com/raerose01/deconstructSigs).

### Run:
```
Rscript -e "rmarkdown::render('analyses/mutational_signatures/run_mutational_signature.Rmd',clean = TRUE)"

```

### Input:
  maf:
  - label: "Input maf file"
  - value: analyses/mutational_signatures/input/pbta-snv-consensus-mutation.maf.tsv.gz
  - input: file

  wgs_bed:
  - label: "Input WGS bed"
  - value: analyses/mutational_signatures/input/WGS.hg38.strelka2.unpadded.bed
  - input: file

  wxs_bed:
  - label: "Input WXS bed"
  - value: analyses/mutational_signatures/input/WXS.hg38.100bp_padded.bed
  - input: file

  metadata_df: 
  - label: "Input metadata"
  - value: analyses/mutational_signatures/input/pbta-histologies.tsv
  - input: file 

  palettes_gradient:
  - label: "Input pallete"
  - value: analyses/mutational_signatures/input/gradient_color_palette.tsv
  - input: file


### Functions:
**run_deconstructSigs()**
Reads in maf file and runs deconstructSigs ( [mut.to.sigs.input()](https://github.com/raerose01/deconstructSigs#muttosigsinput) and  [whichSignatures()](https://github.com/raerose01/deconstructSigs#whichsignatures)) to get a list of weights per signature for each sample.

**bubble_matrix_plot()**
Reads in data.frame from run_deconstructSigs() and plots bubble plots for a given grouping (short_histology [Default])
