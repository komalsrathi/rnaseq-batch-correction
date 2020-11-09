# step1: subset clinical file to RNA-seq and pbta-hgat-dx + create batch column using RNA_library and cohort
Rscript code/02-create-clin.R \
--clin input/pbta-histologies.tsv \
--cohort_col pbta-hgat-dx \
--id_col Kids_First_Biospecimen_ID \
--lib_col RNA_library \
--study_col cohort \
--outfile pbta-hgat-dx_histology_annotation.tsv

# step2: combine matrices
# combine fpkm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds" \
--outfile pbta-gene-expression-rsem-fpkm-collapsed.combined.rds

# combine tpm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-tpm-collapsed.polya.rds" \
--outfile pbta-gene-expression-rsem-tpm-collapsed.combined.rds

# step 3: correct matrices
# correct fpkm matrices
Rscript code/04-batch-correct.R \
--combined_mat output/pbta-gene-expression-rsem-fpkm-collapsed.combined.rds \
--combined_clin input/pbta-hgat-dx_histology_annotation.tsv \
--sample_id Kids_First_Biospecimen_ID \
--corrected_outfile pbta-hgat-dx-gene-expression-rsem-fpkm-corrected.rds \
--uncorrected_outfile pbta-hgat-dx-gene-expression-rsem-fpkm-uncorrected.rds

# correct tpm matrices
Rscript code/04-batch-correct.R \
--combined_mat output/pbta-gene-expression-rsem-tpm-collapsed.combined.rds \
--combined_clin input/pbta-hgat-dx_histology_annotation.tsv \
--sample_id Kids_First_Biospecimen_ID \
--corrected_outfile pbta-hgat-dx-gene-expression-rsem-tpm-corrected.rds \
--uncorrected_outfile pbta-hgat-dx-gene-expression-rsem-tpm-uncorrected.rds

# step4: qc plots 
# corrected tpm matrices
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pbta-hgat-dx-gene-expression-rsem-tpm-uncorrected.rds \
--corrected_mat output/pbta-hgat-dx-gene-expression-rsem-tpm-corrected.rds \
--combined_clin input/pbta-hgat-dx_histology_annotation.tsv \
--sample_id "Kids_First_Biospecimen_ID" \
--tsne_plots pbta-hgat-dx-tpm-tsne.pdf \
--density_plots pbta-hgat-dx-tpm-density.pdf

# corrected fpkm matrices
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pbta-hgat-dx-gene-expression-rsem-fpkm-uncorrected.rds \
--corrected_mat output/pbta-hgat-dx-gene-expression-rsem-fpkm-corrected.rds \
--combined_clin input/pbta-hgat-dx_histology_annotation.tsv \
--sample_id "Kids_First_Biospecimen_ID" \
--tsne_plots pbta-hgat-dx-fpkm-tsne.pdf \
--density_plots pbta-hgat-dx-fpkm-density.pdf


