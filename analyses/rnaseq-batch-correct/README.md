## Batch Correction of TPM Datasets

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to batch correct and combine TPM matrices from various RNA-seq datasets.

Addresses issues: 
1. https://github.com/d3b-center/bixu-tracker/issues/703
2. https://github.com/d3b-center/bixu-tracker/issues/725
3. https://github.com/d3b-center/bixu-tracker/issues/726

### Module Structure

```sh
.
├── README.md
├── code
│   ├── 00-merge-rsem.R			# script to merge RSEM files (if applicable)
│   ├── 01-collapse-matrices.R	# script to uniquify RSEM TPM into gene symbols and sample matrix 
│   ├── 02-combine-matrices.R	# script to combine multiple datasets
│   ├── 02-create-clin.R		# script to create clinical file for each dataset to combine
│   ├── 03-combine-clin.R		# script to combine clinical files created above
│   ├── 04-batch-correct.R		# script to batch correct combined TPM matrices
│   ├── 05-qc-plots.R			# script to create t-SNE and density plots
├── output
│   ├── *-corrected.rds			# matrices for batch corrected datasets
│   ├── *-uncorrected.rds		# matrices for uncorrected datasets
├── plots
│   ├── *-density.pdf			# density plots before and after batch correction
│   └── *-tsne.pdf				# t-SNE plots before and after batch correction
├── run_*_analysis.sh			# script to run full analysis
└── util
    ├── collapse_rnaseq.R		# function to collapse matrices by uniquifying gene symbols
    ├── density_plots.R			# function to create density plots
    ├── install_pkgs.R			# function to install required packages
    ├── pubTheme.R				# function for publication quality ggplot2 theme
    └── tsne_plots.R			# function to create t-SNE plots
```

### Analysis scripts

#### code/00-merge-rsem.R

This script merges all RSEM files into a matrix of gene_id and sample names. Specify --type to create a matrix of FPKM or TPM.

```sh
Rscript code/00-merge-rsem.R --help

Options:
	--rsem_path=RSEM_PATH
		Path to all RSEM genes.results files

	--type=TYPE
		TPM or FPKM

	--outfile=OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/00-merge-rsem.R \
--rsem_path "input/TGEN_normals" \
--type "TPM" \
--outfile "input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds"
```

#### code/01-collapse-matrices.R

This script converts the input matrix (with gene ids) to a matrix of unique gene symbols (row names) and sample identifiers (column names).

```sh
Rscript code/01-collapse-matrices.R --help

Options:
	--mat=MAT
		Expression Matrix (RSEM TPM) (.rds)

	--gene_sym=GENE_SYM
		Is gene symbol present?

	--outfile=OUTFILE
		Output filename (.rds)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/01-collapse-matrices.R \
--mat "input/tgen-brain-normals-gene-expression-rsem-tpm.polya.rds"  \
--gene_sym "FALSE" \
--outfile "input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds"
```

#### code/02-combine-matrices.R

This script will combine all input matrices (with unique gene symbols as row names) into 1 combined matrix. Make sure to run 01-collapse-matrices.R on all individual matrices before combining them.

```sh
Rscript code/02-combine-matrices.R --help

Options:
	--matrices=MATRICES
		Comma separated list of expression matrices to combine (RSEM TPM) (.rds)

	--outfile=OUTFILE
		Output filename (.rds)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/02-combine-matrices.R \
--matrices "input/pbta-gene-expression-rsem-tpm-collapsed.stranded.rds, input/pbta-gene-expression-rsem-tpm-collapsed.polya.rds, input/tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, input/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds" \
--outfile "pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds"
```

#### code/02-create-clin.R

This script will create a clinical or metadata file for each input matrix. If an existing clinical file is used, then it requires the following parameters to be set: `--clin`, `--id_col`, `--lib_col` and `--study_col`. If you don't have any existing clinical file, then use `--mat`, `--lib` and `--study` columns to manually specify the sample ids, library type and study identifier for the dataset.

```sh
Rscript code/02-create-clin.R --help

Options:
	--clin=CLIN
		Existing clinical file (.tsv)

	--cohort_col=COHORT_COL
		cohort column to be used for subsetting clinical file

	--id_col=ID_COL
		Sample identifier column to be used (only use with --clin option)

	--lib_col=LIB_COL
		Library column to be used (only use with --clin option)

	--study_col=STUDY_COL
		Study column to be used (only use with --clin option)

	--mat=MAT
		Expression Matrix (RSEM TPM) (.rds) (when --clin is not provided)

	--lib=LIB
		Library prep method for all samples (when --clin is not provided)

	--study=STUDY
		Study identifier for all samples (when --clin is not provided)

	--outfile=OUTFILE
		Output filename (.tsv)

	-h, --help
		Show this help message and exit

```

##### Example Run

```sh
Rscript code/02-create-clin.R \
--clin "input/pbta-histologies.tsv" \
--cohort_col "pbta-hgat-dx-prog-pm" \
--id_col "Kids_First_Biospecimen_ID" \
--lib_col "RNA_library" \
--study_col "cohort" \
--outfile "pbta-hgat-dx-prog-pm_histology_annotation.tsv"
```

#### code/03-combine-clin.R

This script will combine all input clinical files (generated using code/02-create-clin.R). 

```sh
Rscript code/03-combine-clin.R --help

Options:
	--clin=CLIN
		Comma separated list of metadata to combine (.tsv)

	--outfile=OUTFILE
		Output filename (.tsv)

	-h, --help
		Show this help message and exit
```		

##### Example Run

```sh
Rscript code/03-combine-clin.R \
--clin "input/pbta-hgat-dx-prog-pm_histology_annotation.tsv, input/tgen-batch-metadata.tsv, input/gtex-batch-metadata.tsv" \
--outfile "pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv"
```

#### code/04-batch-correct.R

Script to batch correct combined input dataset generated using code/02-combine-matrices.R. The script using sva::ComBat to batch correct and then back transforms the resulting log2 values for downstream analyses.

##### Example Run

```sh
Rscript code/04-batch-correct.R --help

Options:
	--combined_mat=COMBINED_MAT
		Combined expression matrix with multiple batches (RSEM TPM) (.rds)

	--combined_clin=COMBINED_CLIN
		Combined clinical file with multiple batches (.rds)

	--sample_id=SAMPLE_ID
		Sample identifiers matching between clinical and expression matrix

	--corrected_outfile=CORRECTED_OUTFILE
		Output filename (.RDS)

	--uncorrected_outfile=UNCORRECTED_OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/04-batch-correct.R \
--combined_mat "output/pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds" \
--combined_clin "output/pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv" \
--sample_id "identifier" \
--corrected_outfile "pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-corrected.rds" \
--uncorrected_outfile "pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-uncorrected.rds"
```

#### code/05-qc-plots.R

Script to create QC plots comparing t-SNE and density plots of six house-keeping genes of uncorrected and batch-corrected data.

```sh
Rscript code/05-qc-plots.R --help

Options:
	--uncorrected_mat=UNCORRECTED_MAT
		Combined expression matrix with multiple batches (RSEM TPM) (.rds)

	--corrected_mat=CORRECTED_MAT
		Corrected expression matrix with multiple batches (RSEM TPM) (.rds)

	--combined_clin=COMBINED_CLIN
		Combined clinical file with multiple batches (.tsv)

	--sample_id=SAMPLE_ID
		Sample identifier column in clinical file matching column names in expression datasets

	--tsne_plots=TSNE_PLOTS
		Summary clustering plots (.pdf)

	--density_plots=DENSITY_PLOTS
		Histogram of housekeeping genes (.pdf)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/05-qc-plots.R \
--uncorrected_mat "output/pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-uncorrected.rds" \
--corrected_mat "output/pbta-hgat-dx-prog-pm-tgen-gtex-gene-expression-rsem-tpm-corrected.rds" \
--combined_clin "output/pbta-hgat-dx-prog-pm-tgen-gtex-combined-clin.tsv" \
--sample_id "identifier" \
--tsne_plots "pbta-hgat-dx-prog-pm-tgen-gtex-tpm-tsne.pdf" \
--density_plots "pbta-hgat-dx-prog-pm-tgen-gtex-tpm-density.pdf"
```

### Running the full analysis

```sh
# pbta-hgat-dx-prog-pm + tgen + gtex analysis
bash pbta-hgat-dx-prog-pm-tgen-gtex-analysis.sh
```

#### Output data

All combined uncorrected and corrected matrices can be found here:

```sh
# s3 location:
s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/
```

