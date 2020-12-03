# GTEx Brain
# collapse expected counts 
# collapse tpm
# collapse fpkm  
# these were done separately

# create denovo clinical file
Rscript code/02-create-clin.R \
--mat input/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study GTEx \
--outfile gtex-batch-metadata.tsv