## Coocurrence Plot Generation

Generates plots to display co-occurence and exclusivity of mutations across
multiple tumors and multiple disease types.

### Usage

```
Rscript generate-interaction-plot.R --help

Usage: generate-interaction-plot.R [options]


Options:
	--metadata=CHARACTER
		File path to metadata file to be analyzed. Can be .gz compressed.

	--outdir=OUTDIR
		Output directory path.

	--exclude=EXCLUDE
		File path with a table of genes to be excluded from the figure.
      A tsv file which must contain a column named 'gene` that contains Hugo Symbols

	--samples=CHARACTER
		File path to sample file to be analyzed. Can be .gz compressed.

	--maf=CHARACTER
		File path of MAF file to be analyzed. Can be .gz compressed.

	--palette=CHARACTER
		File path to the palette file for coloring the plots.

	--scripts=SCRIPTS
		Scripts directory path.

	-h, --help
		Show this help message and exit
```

#### Inputs
* A maf file
* A tab separated file that will be the metadata file containing at least the following fields:
	* Kids_First_Biospecimen_ID
	* Kids_First_Participant_ID
	* short_histology
* A tab separated file with fields "Kids_First_Biospecimen_ID" and	"Kids_First_Participant_I"
	* this will be the sample file
* A file specifying the color palette to be used when generating the plots.
* The path to a directory to write the output to.

##### Disease List

The script will analyze all samples together into one co-occurrence plot, and it will also analyze samples and mutations and generate co-ocurrence plots for the following disease list:

* Medulloblastoma
* Low-grade astrocytic tumor
* Ependymoma
* High-grade glioma
* Diffuse midline glioma
* Ganglioglioma
* Craniopharyngioma

#### Outputs

The script outputs a series of plots and files used to generate the plots for each disease type above and for all diseases present in the metadata file. The output directory will have two sub-directories: "figures" which will contain the plots that are generated and "files" which will contain the intermediate files used to generate the plots. Each disease type will have the following files and plots:
* Gene Disease file: tab separated values file with a list of genes and count of samples with mutations in each gene.
* Co-occurrence file: tab separated values file with a list of gene pairs and the co-occurrence score of each pair.
* Co-occurrence plot: plot showing the pairwise co-occurrence of every mutated gene
* Gene disease plot: a bar plot depicted the mutated genes and the number of samples with mutations in each gene.
* Combined plot: both plots combined into a single figure.
