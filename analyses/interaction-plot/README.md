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
      A tsv file which must contain a column named 'gene' that contains Hugo Symbols

	--include=INCLUDE
		File path with a list of genes to be included in the analysis.
      A tsv file which must contain a column named 'gene' that contains Hugo Symbols.
      When using this option, the exclude gene list will be ignored.

	--samples=CHARACTER
		File path to sample file to be analyzed. Can be .gz compressed.

	--maf=CHARACTER
		File path of MAF file to be analyzed. Can be .gz compressed.

	--palette=CHARACTER
		File path to the palette file for coloring the plots.

	--scripts=SCRIPTS
		Scripts directory path.

	--p_cut_off=P_CUT_OFF
		Highest allowable p value for cooccurrence pairs. Default value: 1,
      no filtering

	--q_cut_off=Q_CUT_OFF
		q value cut off for cooccurrence pairs. Default value: 1,
      no filtering

	--include_syn
		Include synonymous coding mutations

	--include_noncoding
		Include noncoding mutations (within transcript)

	--include_nontranscribed
		Include nontranscribed (upstream & downstream) mutations

	-h, --help
		Show this help message and exit
```

#### Inputs
* A maf file
* A tab separated file that will be the metadata file containing at least the following fields:
	* Kids_First_Biospecimen_ID
	* Kids_First_Participant_ID
	* short_histology
* A tab separated file with fields "Kids_First_Biospecimen_ID" and	"Kids_First_Participant_ID"
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

#### Example Command

```
Rscript generate-interaction-plot.R --palette divergent_color_palette.tsv --q_cut_off 0.05 --p_cut_off 0.05 --include ./gene_lists/all_gene_list.tsv --metadata /data/pnoc003_histologies_v17_candidate.tsv --samples /data/pnoc003_wgs_samples.tsv --maf /data/pnoc003_wgs-consensus-open-pbta-merged.maf.gz --out /data/plot-out/p05-inc/pnoc003-wgs
```

#### Outputs

The script outputs a series of plots and files used to generate the plots for each disease type above and for all diseases present in the metadata file. The output directory will have two sub-directories: "figures" which will contain the plots that are generated and "files" which will contain the intermediate files used to generate the plots. Each disease type will have the following files and plots:
* Gene Disease file: tab separated values file with a list of genes and count of samples with mutations in each gene.
* Co-occurrence file: tab separated values file with a list of gene pairs and the co-occurrence score of each pair.
* Co-occurrence plot: plot showing the pairwise co-occurrence of every mutated gene
* Correlation plot: Co-occurrence plot: plot showing the pairwise co-occurrence of every mutated gene. This may replace the co-occurrence plot in future updates.
* Gene disease plot: a bar plot depicted the mutated genes and the number of samples with mutations in each gene.
* Combined plot: both plots combined into a single figure.

#### Key Functions

The main analysis is carried out by helper functions found in the scripts subdirectory and controlled from the generate-interaction-plot.R script.

The key functions are:

* gene_counts() takes a maf file with mutations from samples being analyzed and uses dplyr functions to generate a count of genes with mutations for each sample in the maf file.
* coocurrence() calculates the cooccurrence or mutual exclusivity of gene pairs and assigns each pair a significance value based on Fisher's exact test.
* modify_cooc_sum() removes insignificant cooccurrence pairs based off a provided p value
* plot_corr() uses the R corrplot library to graph the cooccurrence scores and label pairs with a hypothesis adjused value greater than a provided value
* plot_disease() creates a bar plot with a count of mutations in each gene and disease type
