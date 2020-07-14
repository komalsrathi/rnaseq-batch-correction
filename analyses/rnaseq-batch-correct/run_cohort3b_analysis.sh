# 01. Collapse RNA-seq
# PNOC003 Cohort3b polya
Rscript code/01-collapse-matrices.R \
--mat input/cohort3b_subset-gene-expression-rsem-tpm.polya.rds \
--gene_sym FALSE \
--outfile cohort3b_subset-gene-expression-rsem-tpm.collapsed.polya.rds

# PNOC003 Cohort3b stranded
Rscript code/01-collapse-matrices.R \
--mat input/cohort3b_subset-gene-expression-rsem-tpm.stranded.rds \
--gene_sym FALSE \
--outfile cohort3b_subset-gene-expression-rsem-tpm.collapsed.stranded.rds

# 02. Combine collapsed matrices (Cohort 3b + GTEx + TGEN)
Rscript code/02-combine-matrices.R \
--matrices "output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.polya.rds, output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.stranded.rds, output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, output/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds" \
--outfile pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm.rds

# 02. Combine collapsed matrices for immune deconvolution (Cohort 3b polyA + stranded)
Rscript code/02-combine-matrices.R \
--matrices "output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.polya.rds, output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.stranded.rds" \
--outfile pnoc003-cohort3b-gene-expression-rsem-tpm.rds

# 02. Create Clinical for each study
# PNOC003 Cohort3b polyA
Rscript code/02-create-clin.R \
--mat output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.polya.rds \
--clin input/hgat_all_primary_longitudinal.rds \
--id_col 'Kids_First_Biospecimen_ID' \
--lib_col 'library' \
--study_col 'cohort' \
--outfile pnoc003-cohort3b-polya-batch-metadata.rds

# PNOC003 Cohort3b stranded
Rscript code/02-create-clin.R \
--mat output/cohort3b_subset-gene-expression-rsem-tpm.collapsed.stranded.rds \
--clin input/hgat_all_primary_longitudinal.rds \
--id_col 'Kids_First_Biospecimen_ID' \
--lib_col 'library' \
--study_col 'cohort' \
--outfile pnoc003-cohort3b-stranded-batch-metadata.rds

# 03. Combine Clinical (Cohort 3b + GTEx + TGEN)
Rscript code/03-combine-clin.R \
--clin 'input/pnoc003-cohort3b-polya-batch-metadata.rds, input/pnoc003-cohort3b-stranded-batch-metadata.rds, input/gtex-batch-metadata.rds, input/tgen-batch-metadata.rds' \
--outfile pnoc003-cohort3b-gtex-tgen-metadata.rds

# 04. Batch Correction (Cohort 3b + GTEx + TGEN)
Rscript code/04-batch-correct.R \
--combined_mat output/pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm.rds \
--combined_clin output/pnoc003-cohort3b-gtex-tgen-metadata.rds \
--outfile pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# 05. QC Plots (Cohort 3b + GTEx + TGEN)
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm.rds \
--corrected_mat output/pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm-corrected.rds \
--combined_clin output/pnoc003-cohort3b-gtex-tgen-metadata.rds \
--tsne_plots pnoc003-cohort3b-gtex-tgen-tsne.pdf \
--density_plots pnoc003-cohort3b-gtex-tgen-density.pdf
