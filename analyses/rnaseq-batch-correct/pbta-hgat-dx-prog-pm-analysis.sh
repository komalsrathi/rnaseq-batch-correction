# step1: subset clinical file to RNA-seq and pbta-hgat-dx-prog-pm + create batch column using RNA_library and cohort
# create clinical file for tpm based batch correction
# for tpm based batch correction, we are using batch as a combination of study_id + rna_library
Rscript code/02-create-clin.R \
--clin "input/pbta-histologies.tsv" \
--cohort_col "pbta-hgat-dx-prog-pm" \
--id_col "Kids_First_Biospecimen_ID" \
--lib_col "RNA_library" \
--study_col "cohort" \
--outfile "pbta-hgat-dx-prog-pm_histology_annotation.tsv"

# adding the two additional samples manually: PNOC008-23 (BS_00FD2KMP) and PT_SMJNYRBK (BS_ASM7ANKQ)
# BS_00FD2KMP	PNOC008_rna_exome
# BS_ASM7ANKQ	CBTN_stranded
echo "BS_00FD2KMP	PNOC008_rna_exome" >> input/pbta-hgat-dx-prog-pm_histology_annotation.tsv
echo "BS_ASM7ANKQ	CBTN_stranded" >> input/pbta-hgat-dx-prog-pm_histology_annotation.tsv

# create clinical file for count based batch correction
# for count based batch correction, we are using batch only with rna_library. Issue discussed here: https://github.com/d3b-center/bixu-tracker/issues/875#issuecomment-738814246
Rscript code/02-create-clin.R \
--clin "input/pbta-histologies.tsv" \
--cohort_col "pbta-hgat-dx-prog-pm" \
--id_col "Kids_First_Biospecimen_ID" \
--lib_col "RNA_library" \
--study_col "" \
--outfile "pbta-hgat-dx-prog-pm_histology_annotation-forcounts.tsv"

# adding the two additional samples manually: PNOC008-23 (BS_00FD2KMP) and PT_SMJNYRBK (BS_ASM7ANKQ)
# BS_00FD2KMP	rna_exome
# BS_ASM7ANKQ	stranded
echo "BS_00FD2KMP	rna_exome" >> input/pbta-hgat-dx-prog-pm_histology_annotation-forcounts.tsv
echo "BS_ASM7ANKQ	stranded" >> input/pbta-hgat-dx-prog-pm_histology_annotation-forcounts.tsv

# step2: combine matrices
# combine tpm matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-tpm-collapsed.polya.rds" \
--outfile "pbta-gene-expression-rsem-tpm-collapsed.combined.rds"

# combine expected count matrices
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds, input/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds" \
--outfile "pbta-gene-counts-rsem-expected_count-collapsed.combined.rds"

# step3: combine clinical
# here we can use the same clinical file as generated in step1 (no need to combine any clinical files)

# step 4: batch-correct matrices
# correct tpm matrices
Rscript code/04-batch-correct.R \
--combined_mat "output/pbta-gene-expression-rsem-tpm-collapsed.combined.rds" \
--combined_clin "input/pbta-hgat-dx-prog-pm_histology_annotation.tsv" \
--sample_id "identifier" \
--type "TPM" \
--corrected_outfile "pbta-hgat-dx-prog-pm-gene-expression-rsem-tpm-corrected.rds" \
--uncorrected_outfile "pbta-hgat-dx-prog-pm-gene-expression-rsem-tpm-uncorrected.rds"

# correct count matrices
Rscript code/04-batch-correct.R \
--combined_mat "output/pbta-gene-counts-rsem-expected_count-collapsed.combined.rds" \
--combined_clin "input/pbta-hgat-dx-prog-pm_histology_annotation-forcounts.tsv" \
--sample_id "identifier" \
--type "expected_count" \
--corrected_outfile "pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-corrected.rds" \
--uncorrected_outfile "pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds"

# step 5: qc plots
# corrected tpm matrices
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pbta-hgat-dx-prog-pm-gene-expression-rsem-tpm-uncorrected.rds \
--corrected_mat output/pbta-hgat-dx-prog-pm-gene-expression-rsem-tpm-corrected.rds \
--combined_clin input/pbta-hgat-dx-prog-pm_histology_annotation.tsv \
--sample_id "identifier" \
--tsne_plots pbta-hgat-dx-prog-pm-tpm-tsne.pdf \
--density_plots pbta-hgat-dx-prog-pm-tpm-density.pdf

# corrected count matrices
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds \
--corrected_mat output/pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-corrected.rds \
--combined_clin input/pbta-hgat-dx-prog-pm_histology_annotation-forcounts.tsv \
--sample_id "identifier" \
--tsne_plots pbta-hgat-dx-prog-pm-expected_count-tsne.pdf \
--density_plots pbta-hgat-dx-prog-pm-expected_count-density.pdf
