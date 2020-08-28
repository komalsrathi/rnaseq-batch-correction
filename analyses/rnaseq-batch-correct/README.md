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
│   ├── 00-data-prep.R			# script to merge RSEM files (if applicable)
│   ├── 01-collapse-matrices.R	# script to uniquify RSEM TPM into gene symbols and sample matrix 
│   ├── 02-combine-matrices.R	# script to combine multiple datasets
│   ├── 02-create-clin.R		# script to create clinical file for each dataset to combine
│   ├── 03-combine-clin.R		# script to combine clinical files created above
│   ├── 04-batch-correct.R		# script to batch correct combined TPM matrices
│   ├── 05-qc-plots.R			# script to create t-SNE and density plots
├── output
│   ├── *-corrected.rds			# combined matrices for batch corrected datasets
	├── *-metadata.rds			# combined clinical file for batch corrected datasets
├── plots
│   ├── *-density.pdf			# density plots before and after batch correction
│   └── *-tsne.pdf				# t-SNE plots before and after batch correction
├── run_*_analysis.sh			# script to run full analysis
└── util
    ├── collapse_rnaseq.R		# function to collapse matrices by uniquifying gene symbols
    ├── density_plots.R			# function to create density plots
    ├── pubTheme.R				# function for publication quality ggplot2 theme
    └── tsne_plots.R			# function to create t-SNE plots
```

### Analysis scripts

#### code/00-data-prep.R

This script merges all RSEM files into a matrix of gene_id and sample names.

```sh
Rscript code/00-data-prep.R --help

Options:
	--rsem_path=RSEM_PATH
		Path to all RSEM genes.results files

	--outfile=OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/00-data-prep.R \
--rsem_path input/TGEN_normals \
--outfile tgen-normals-gene-expression-rsem-tpm.polya.rds
```

#### code/01-collapse-matrices.R

This script converts the input matrix (with gene ids) to a matrix of unique gene symbols (row names) and sample identifiers (column names).

```sh
Rscript code/01-collapse-matrices.R --help

Options:
	--mat=MAT
		Expression Matrix (RSEM TPM) (.RDS)

	--gene_sym=GENE_SYM
		Is gene symbol present?

	--outfile=OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/01-collapse-matrices.R \
--mat input/tgen-normals-gene-expression-rsem-tpm.polya.rds \
--gene_sym FALSE \
--outfile tgen-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds
```

#### code/02-combine-matrices.R

This script will combine all input matrices (with unique gene symbols as row names) into 1 combined matrix. Make sure to run 01-collapse-matrices.R on all individual matrices before combining them.

```sh
Rscript code/02-combine-matrices.R --help

Options:
	--matrices=MATRICES
		Comma separated list of expression matrices to combine (RSEM TPM) (.RDS)

	--outfile=OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
Rscript code/02-combine-matrices.R \
--matrices "output/pnoc003_subset_Diagnosis-gene-expression-rsem-tpm.collapsed.polya.rds, output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds, output
--outfile pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds
```

#### code/02-create-clin.R

This script will create a clinical or metadata file for each input matrix. If an existing clinical file is used, then it requires the following parameters to be set: `--clin`, `--id_col`, `--lib_col` and `--study_col`. If you don't have any existing clinical file, then use `--lib` and `--study` columns to manually specify the library type and study identifier for the dataset.

```sh
Rscript code/02-create-clin.R --help

Options:
	--mat=MAT
		Expression Matrix (RSEM TPM) (.rds)

	--clin=CLIN
		Existing clinical file (.rds)

	--id_col=ID_COL
		Sample identifier column to be used (only use with --clin option)

	--lib_col=LIB_COL
		Library column to be used (only use with --clin option)

	--study_col=STUDY_COL
		Study column to be used (only use with --clin option)

	--lib=LIB
		Library Prep Method

	--study=STUDY
		Study identifier

	--outfile=OUTFILE
		Output filename (.rds)

	-h, --help
		Show this help message and exit
```

##### Example Run

```sh
# when clinical file is available
Rscript code/02-create-clin.R \
--mat output/cohort3a_subset-gene-expression-rsem-tpm.collapsed.polya.rds \
--clin input/hgat_all_primary.rds \
--id_col 'Kids_First_Biospecimen_ID' \
--lib_col 'library' \
--study_col 'cohort' \
--outfile pnoc003-cohort3a-polya-batch-metadata.rds

# when clinical file is not available
Rscript code/02-create-clin.R \
--mat output/gtex-brain-normals-gene-expression-rsem-tpm.collapsed.polya.rds \
--lib polyA \
--study GTEx \
--outfile gtex-batch-metadata.rds
```

#### code/03-combine-clin.R

This script will combine all input clinical files (generated using code/02-create-clin.R). 

```sh
Rscript code/03-combine-clin.R --help

Options:
	--clin=CLIN
		Comma separated list of metadata to combine (.RDS)

	--outfile=OUTFILE
		Output filename (.RDS)

	-h, --help
		Show this help message and exit
```		

##### Example Run

```sh
Rscript code/03-combine-clin.R \
--clin 'input/pnoc003-cohort1-batch-metadata.rds, input/gtex-batch-metadata.rds, input/tgen-batch-metadata.rds' \
--outfile pnoc003-cohort1-gtex-tgen-metadata.rds
```

#### code/04-batch-correct.R

Script to batch correct combined input dataset generated using code/02-combine-matrices.R. The script using sva::ComBat to batch correct and then back transforms the resulting log2 values for downstream analyses.

##### Example Run

```sh
Rscript code/04-batch-correct.R --help

Rscript code/04-batch-correct.R \
--combined_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds \
--combined_clin output/pnoc003-cohort1-gtex-tgen-metadata.rds \
--outfile pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm-corrected.rds
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
		Combined clinical file with multiple batches (.rds)

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
--uncorrected_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds \
--corrected_mat output/pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm-corrected.rds \
--combined_clin output/pnoc003-cohort1-gtex-tgen-metadata.rds \
--tsne_plots pnoc003-cohort1-gtex-tgen-tsne.pdf \
--density_plots pnoc003-cohort1-gtex-tgen-density.pdf
```

### Running the full analysis

```sh
# Use any of the following scripts as a template to run the full analysis
# PNOC003 Cohort 1 + GTEx Brain + TGEN Brain
bash run_cohort1_analysis.sh

# PNOC003 Cohort 3a + GTEx Brain + TGEN Brain
bash run_cohort3a_analysis.sh

# PNOC003 Cohort 3b + GTEx Brain + TGEN Brain
bash run_cohort3b_analysis.sh
```

#### Output data

All combined uncorrected and corrected matrices, combined clinical files and qc plots can be found here:

```sh
# s3 location:
s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/batch_correction/

# command used
aws s3 --profile saml sync output/ s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/batch_correction/ --include "*.rds"
aws s3 --profile saml sync plots/ s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/batch_correction/ --include "*.pdf"
```

Example uncorrected matrices:

```sh
# PNOC003 Cohort 1 + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm.rds

# PNOC003 Cohort 3a + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort3a-gtex-tgen-gene-expression-rsem-tpm.rds

# PNOC003 Cohort 3a (no PNOC008) + GTEx Normal Brain + TGEN Brain :
pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm.rds

# PNOC003 Cohort 3b + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm.rds
```

Example corrected matrices:

```sh
# PNOC003 Cohort 1 + GTEx Normal Brain + TGEN Brain: 
pnoc003-cohort1-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# PNOC003 Cohort 3a + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort3a-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# PNOC003 Cohort 3a (no PNOC008) + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort3a-wo-PNOC008-gtex-tgen-gene-expression-rsem-tpm-corrected.rds

# PNOC003 Cohort 3b + GTEx Normal Brain + TGEN Brain : 
pnoc003-cohort3b-gtex-tgen-gene-expression-rsem-tpm-corrected.rds
```