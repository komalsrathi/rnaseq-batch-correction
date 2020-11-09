# step1: subset clinical file to RNA-seq and pbta-hgat-dx-prog-pm + create batch column using RNA_library and cohort
Rscript code/02-create-clin.R \
--clin "input/pbta-histologies.tsv" \
--cohort_col "pbta-hgat-dx-prog-pm" \
--id_col "Kids_First_Biospecimen_ID" \
--lib_col "RNA_library" \
--study_col "cohort" \
--outfile "pbta-hgat-dx-prog-pm_histology_annotation.tsv"

# step2: combine matrices
# combine tpm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-tpm-collapsed.polya.rds, input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, input/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds" \
--outfile "pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds"

# step3: combine clinical
Rscript code/03-combine-clin.R \
--clin "input/pbta-hgat-dx-prog-pm_histology_annotation.tsv, input/tgen-batch-metadata.tsv, input/gtex-batch-metadata.tsv" \
--outfile "pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv"

# step 4: correct matrices
# correct tpm matrices
Rscript code/04-batch-correct.R \
--combined_mat "output/pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds" \
--combined_clin "output/pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv" \
--sample_id "identifier" \
--corrected_outfile "pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-corrected.rds" \
--uncorrected_outfile "pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-uncorrected.rds"

# step 5: qc plots
# corrected tpm matrices
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-uncorrected.rds \
--corrected_mat output/pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-corrected.rds \
--combined_clin output/pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv \
--sample_id "identifier" \
--tsne_plots pbta-hgat-dx-prog-pm-tgen-gtex-tpm-tsne.pdf \
--density_plots pbta-hgat-dx-prog-pm-tgen-gtex-tpm-density.pdf