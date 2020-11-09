# GTEx Brain
# collapse RNA-seq TPM
Rscript code/01-collapse-matrices.R \
--mat input/gtex-brain-normals-gene-expression-rsem-tpm.polya.rds  \
--gene_sym TRUE \
--outfile input/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds

# create denovo clinical file
Rscript code/02-create-clin.R \
--mat input/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study GTEx \
--outfile gtex-batch-metadata.tsv