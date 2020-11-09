# TGEN normals
# merge RSEM TPM
Rscript code/00-merge-rsem.R \
--rsem_path input/TGEN_normals \
--type TPM \
--outfile input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds

# collapse RNA-seq TPM
Rscript code/01-collapse-matrices.R \
--mat input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds  \
--gene_sym FALSE \
--outfile input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds

# create denovo clinical
Rscript code/02-create-clin.R \
--mat input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study TGEN \
--outfile tgen-batch-metadata.tsv