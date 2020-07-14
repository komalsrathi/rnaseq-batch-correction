# 00. Data prep
Rscript code/00-data-prep.R \
--rsem_path input/TGEN_normals \
--outfile tgen-normals-gene-expression-rsem-tpm.polya.rds

# 01. Collapse RNA-seq
# GTEx Brain
Rscript code/01-collapse-matrices.R \
--mat input/gtex-brain-normals-gene-expression-rsem-tpm.polya.rds  \
--gene_sym TRUE \
--outfile gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds

# TGEN Brain
Rscript code/01-collapse-matrices.R \
--mat input/tgen-normals-gene-expression-rsem-tpm.polya.rds \
--gene_sym FALSE \
--outfile tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds

# PNOC003 Cohort1
Rscript code/01-collapse-matrices.R \
--mat input/pnoc003_subset_Diagnosis-gene-expression-rsem-tpm.polya.rds \
--gene_sym FALSE \
--outfile pnoc003_subset_Diagnosis-gene-expression-rsem-tpm.collapsed.polya.rds

# 02. Combine collapsed matrices
Rscript code/02-combine-matrices.R \
--matrices "output/pnoc003_subset_Diagnosis-gene-expression-rsem-tpm.collapsed.polya.rds, output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, output/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds" \
--outfile pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds

# 02. Create Clinical for each study
# GTEx Brain
Rscript code/02-create-clin.R \
--mat output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study GTEx \
--outfile gtex-batch-metadata.rds

# TGEN Brain
Rscript code/02-create-clin.R \
--mat output/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study TGEN \
--outfile tgen-batch-metadata.rds

# PNOC003 Cohort1
Rscript code/02-create-clin.R \
--mat output/pnoc003_subset_Diagnosis-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study PNOC003_cohort1 \
--outfile pnoc003-cohort1-batch-metadata.rds

# 03. Combine Clinical
Rscript code/03-combine-clin.R \
--clin 'input/pnoc003-cohort1-batch-metadata.rds, input/gtex-batch-metadata.rds, input/tgen-batch-metadata.rds' \
--outfile pnoc003-cohort1-gtex-tgen-metadata.rds

# 04. Batch Correction
Rscript code/04-batch-correct.R \
--combined_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds \
--combined_clin output/pnoc003-cohort1-gtex-tgen-metadata.rds \
--outfile pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# 05. QC Plots
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds \
--corrected_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm-corrected.rds \
--combined_clin output/pnoc003-cohort1-gtex-tgen-metadata.rds \
--tsne_plots pnoc003-cohort1-gtex-tgen-tsne.pdf \
--density_plots pnoc003-cohort1-gtex-tgen-density.pdf
