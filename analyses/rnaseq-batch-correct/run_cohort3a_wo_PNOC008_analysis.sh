# 01. Collapse RNA-seq
# PNOC003 Cohort3a polya
echo "Collapse RNA-seq polya"
Rscript code/01-collapse-matrices.R \
--mat input/cohort3a_subset-gene-expression-rsem-tpm.polya.rds \
--gene_sym FALSE \
--outfile cohort3a_subset-gene-expression-rsem-tpm.collapsed.polya.rds

# PNOC003 Cohort3a stranded (w/o PNOC008)
echo "Collapse RNA-seq stranded"
Rscript code/01-collapse-matrices.R \
--mat input/cohort3a_subset-wo-PNOC008-gene-expression-rsem-tpm.stranded.rds \
--gene_sym FALSE \
--outfile cohort3a_subset-wo-PNOC008-gene-expression-rsem-tpm.collapsed.stranded.rds

# 02. Combine collapsed matrices (Cohort 3a w/o PNOC008 + GTEx + TGEN)
echo "Combine collapsed matrices"
Rscript code/02-combine-matrices.R \
--matrices "output/cohort3a_subset-gene-expression-rsem-tpm.collapsed.polya.rds, output/cohort3a_subset-wo-PNOC008-gene-expression-rsem-tpm.collapsed.stranded.rds, output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, output/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds" \
--outfile pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm.rds

# 02. Combine collapsed matrices for immune deconvolution (Cohort 3a polyA + stranded)
echo "Combine collapsed matrices (just 3a)"
Rscript code/02-combine-matrices.R \
--matrices "output/cohort3a_subset-gene-expression-rsem-tpm.collapsed.polya.rds, output/cohort3a_subset-wo-PNOC008-gene-expression-rsem-tpm.collapsed.stranded.rds" \
--outfile pnoc003-cohort3a-wo-PNOC008-gene-expression-rsem-tpm.rds

# 02. Create Clinical for each study
# PNOC003 Cohort3a polyA
echo "Create clinical file: polyA"
Rscript code/02-create-clin.R \
--mat output/cohort3a_subset-gene-expression-rsem-tpm.collapsed.polya.rds \
--clin input/hgat_all_primary.rds \
--id_col 'Kids_First_Biospecimen_ID' \
--lib_col 'library' \
--study_col 'cohort' \
--outfile pnoc003-cohort3a-polya-batch-metadata.rds

# PNOC003 Cohort3a stranded
echo "Create clinical file: stranded"
Rscript code/02-create-clin.R \
--mat output/cohort3a_subset-wo-PNOC008-gene-expression-rsem-tpm.collapsed.stranded.rds \
--clin input/hgat_all_primary.rds \
--id_col 'Kids_First_Biospecimen_ID' \
--lib_col 'library' \
--study_col 'cohort' \
--outfile pnoc003-cohort3a-wo-PNOC008-stranded-batch-metadata.rds

# 03. Combine Clinical (Cohort 3a + GTEx + TGEN)
echo "Create combined clinical file"
Rscript code/03-combine-clin.R \
--clin 'input/pnoc003-cohort3a-polya-batch-metadata.rds, input/pnoc003-cohort3a-wo-PNOC008-stranded-batch-metadata.rds, input/gtex-batch-metadata.rds, input/tgen-batch-metadata.rds' \
--outfile pnoc003-cohort3a-wo-PNOC008-gtex-tgen-metadata.rds

# 04. Batch Correction (Cohort 3a + GTEx + TGEN)
echo "Batch correction"
Rscript code/04-batch-correct.R \
--combined_mat output/pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm.rds \
--combined_clin output/pnoc003-cohort3a-wo-PNOC008-gtex-tgen-metadata.rds \
--outfile pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# 05. QC Plots (Cohort 3a + GTEx + TGEN)
echo "QC Plots"
Rscript code/05-qc-plots.R \
--uncorrected_mat output/pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm.rds \
--corrected_mat output/pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm-corrected.rds \
--combined_clin output/pnoc003-cohort3a-wo-PNOC008-gtex-tgen-metadata.rds \
--tsne_plots pnoc003-cohort3a-wo-PNOC008-gtex-tgen-tsne.pdf \
--density_plots pnoc003-cohort3a-wo-PNOC008-gtex-tgen-density.pdf
