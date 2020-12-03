# TGEN normals

# merge rsem files
# tpm matrix
Rscript code/00-merge-rsem.R \
--rsem_path input/TGEN_normals \
--type TPM \
--outfile input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds

# fpkm matrix
Rscript code/00-merge-rsem.R \
--rsem_path input/TGEN_normals \
--type FPKM \
--outfile input/tgen-brain-normals-gene-expression-rsem-fpkm.polya.rds

# expected count matrix
Rscript code/00-merge-rsem.R \
--rsem_path input/TGEN_normals \
--type expected_count \
--outfile input/tgen-brain-normals-gene-counts-rsem-expected_count.polya.rds

# collapse to gene symbols
# collapse tpm
Rscript code/01-collapse-matrices.R \
--mat input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds  \
--gene_sym FALSE \
--outfile input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds

# collapse fpkm
Rscript code/01-collapse-matrices.R \
--mat input/tgen-brain-normals-gene-expression-rsem-fpkm.polya.rds  \
--gene_sym FALSE \
--outfile input/tgen-brain-normals-gene-expression-rsem-fpkm.collapsed.polya.rds

# collapse expected counts
Rscript code/01-collapse-matrices.R \
--mat input/tgen-brain-normals-gene-counts-rsem-expected_count.polya.rds  \
--gene_sym FALSE \
--outfile input/tgen-brain-normals-gene-counts-rsem-expected_count.collapsed.polya.rds

# create denovo clinical
Rscript code/02-create-clin.R \
--mat input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study TGEN \
--outfile tgen-batch-metadata.tsv