## Coocurrence Plot Generation

Generates plots to display co-occurence and exclusivity of mutations across
multiple tumors and multiple disease types.

### Usage

```
Rscript generate-interaction-plot.R --help

Usage: generate-interaction-plot.R [options]


Options:
	--metadata=CHARACTER
		File path to MAF file to be analyzed. Can be .gz compressed.

	--outdir=OUTDIR
		Output directory path.

	--exclude=EXCLUDE
		File path with a table of genes to be excluded from the figure.
      A tsv file which must contain a column named 'gene` that contains Hugo Symbols

	--samples=CHARACTER
		File path to MAF file to be analyzed. Can be .gz compressed.

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

#### Outputs
