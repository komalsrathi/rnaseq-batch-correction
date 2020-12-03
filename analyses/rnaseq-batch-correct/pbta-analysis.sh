# combine pbta matrices
# combine fpkm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds" \
--outfile pbta-gene-expression-rsem-fpkm-collapsed.combined.rds

# combine tpm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-tpm-collapsed.polya.rds" \
--outfile pbta-gene-expression-rsem-tpm-collapsed.combined.rds

# combine expected count matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds, input/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds" \
--outfile pbta-gene-counts-rsem-expected_count-collapsed.combined.rds