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

### Functions:
**run_deconstructSigs()**
Reads in maf file and runs deconstructSigs ( [mut.to.sigs.input()](https://github.com/raerose01/deconstructSigs#muttosigsinput) and  [whichSignatures()](https://github.com/raerose01/deconstructSigs#whichsignatures)) to get a list of weights per signature for each sample.

