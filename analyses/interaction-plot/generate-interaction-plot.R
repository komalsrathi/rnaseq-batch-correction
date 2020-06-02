#R script wrapper for generating interaction plots

#TODO / Questions:
# Add CNV options and parsing?

#Magrittr pipe
#symbol set to basically work like unix pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)
library(ggplot2)
library(hash)
library(patchwork)

option_list <- list(
  make_option(
    opt_str = "--metadata",
    type = "character",
    help = "File path to metadata file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outdir",
    default = file.path(getwd(), "out"),
    type = "character",
    help = "Output directory path."
  ),
  make_option(
    opt_str = "--exclude",
    default = file.path(getwd(), "exclude-genes.txt"),
    type = "character",
    help = "File path with a table of genes to be excluded from the figure.
      A tsv file which must contain a column named 'gene' that contains Hugo Symbols"
  ),
  make_option(
    opt_str = "--include",
    type = "character",
    help = "File path with a list of genes to be included in the analysis.
      A tsv file which must contain a column named 'gene' that contains Hugo Symbols.
      When using this option, the exclude gene list will be ignored."
  ),
  make_option(
    opt_str = "--samples",
    type = "character",
    help = "File path to sample file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--maf",
    type = "character",
    help = "File path of MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--palette",
    type = "character",
    help = "File path to the palette file for coloring the plots.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--scripts",
    default = file.path(getwd(), "scripts"),
    type = "character",
    help = "Scripts directory path."
  ),
  make_option(
    opt_str = "--p_cut_off",
    default = 1,
    type = "numeric",
    help = "Highest allowable p value for cooccurrence pairs. Default value: 1,
      no filtering"
  ),
  make_option(
    opt_str = "--q_cut_off",
    default = 1,
    type = "numeric",
    help = "q value cut off for cooccurrence pairs. Default value: 1,
      no filtering"
  ),
  make_option(
  opt_str = "--include_syn",
  action = "store_true",
  default = FALSE,

  help = "Include synonymous coding mutations"
  ),
  make_option(
    opt_str = "--include_noncoding",
    action = "store_true",
    default = FALSE,
    help = "Include noncoding mutations (within transcript)"
  ),
  make_option(
    opt_str = "--include_nontranscribed",
    action = "store_true",
    default = FALSE,
    help = "Include nontranscribed (upstream & downstream) mutations"
  )
)

#defaults that might be changed later
max_genes <- 50
min_mutated <- 5 #Number of samples that have to have mutation for it to be included
plot_size <- 50

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
meta_file <- file.path(opts$metadata)
sample_file <- file.path(opts$samples)
maf_file <- file.path(opts$maf)
out_dir <- file.path(opts$outdir)
exclude_file <- file.path(opts$exclude)
script_dir <- file.path(opts$scripts)
palette_file <- file.path(opts$palette)
p_cut <- opts$p_cut_off
q_cut <- opts$q_cut_off
print(meta_file)
print(out_dir)
print(script_dir)

#load additional / helper functions
source(file.path(script_dir, "cooccur_functions.R"))
source(file.path(script_dir, "process_inputs.R"))
source(file.path(script_dir, "reduce-maf.R"))
source(file.path(script_dir, "calculate_score.R"))
source(file.path(script_dir, "make_plots.R"))

#set up color palette
colors <- set_colors(palette_file)
divergent_colors <- colors$divergent_colors
na_color <- colors$na_color

#create output directories
figure_dir <- file.path(out_dir, "figures")
file_dir <- file.path(out_dir, "files")
dir.create(out_dir)
dir.create(figure_dir)
dir.create(file_dir)

#set up diseases
diseases <- hash(
"all" = "All",
"Medulloblastoma" = "Medulloblastoma",
"LGAT" = "Low-grade astrocytic tumor",
"Ependymoma" = "Ependymoma",
"HGAT" = "High-grade glioma",
"DMG" = "Diffuse midline glioma",
"Ganglioglioma" = "Ganglioglioma",
"Craniopharyngioma" = "Craniopharyngioma"
)

#determine which variant types to include
include_types <- set_var_types(opts)

#read inputs
meta_df <- readr::read_tsv(meta_file, col_types = readr::cols(),
  guess_max = 10000) %>%
  dplyr::select("Kids_First_Biospecimen_ID","short_histology") %>%
  unique()

sample_df <- readr::read_tsv(sample_file, col_types = readr::cols())
maf_df <- data.table::fread(maf_file, data.table = FALSE)
maf_df <- reduce_maf(maf_df)

if (! is.null(opts$include)) {
  include_file <- file.path(opts$include)
  include_list <- readr::read_tsv(include_file, col_types = readr::cols())
}

#run analysis for each disease
for (disease_id in keys(diseases)) {
  print(disease_id)

  #get the gene list
  if (exists("include_list")) {
    genes <- include_list$gene
  } else {
    genes <- get_gene_list(maf_df, exclude_file, disease_id)
  }

  #get samples with the disease
  samples <- samples_with_disease(meta_df, sample_df, disease_id)

  #reduce meta_df to only the needed samples
  sample_meta <- meta_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% samples)

  #filter maf to chosen samples
  maf_filtered <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% samples)

  #don't run if there's no samples to run
  row_count <- nrow(maf_filtered)
  if (row_count == 0) {
    print(paste("No samples to be analyzed for", disease_id))
    next
  }

  #remove unwanted variant types from maf
  maf_filtered <- filter_mutations(maf_filtered, include_types)

  # count mutations by gene/sample pair
  gene_sample_counts <- gene_counts(maf_filtered, genes)

  #get most often mutated genes
  top_count_genes <- get_top_genes(gene_sample_counts, max_genes, min_mutated)

  #count mutated samples by disease type
  gene_disease_counts <- count_diseases(gene_sample_counts, top_count_genes,
    sample_meta)

  #calculate coocurrence
  cooccur_summary <- coocurrence(gene_sample_counts, top_count_genes)

  #remove insignificant gene pairs from cooccur_summary
  cooccur_sig <- sig_filter(cooccur_summary, p_cut)

  #don't run if there's no significant coccurrence pairs to run
  row_count <- nrow(cooccur_sig)
  if (row_count == 0) {
    print(paste("No significant cooccurrences for", disease_id))
    next
  }

  #reduce cooccur summary to needed fields and generate labels
  coocur_df <- modify_cooc_sum(cooccur_sig)

  #try out corrplot
  corr_file <- file.path(figure_dir, paste("corr.", disease_id, ".png",
    sep = ""))
  corr_plot <- plot_corr(coocur_df, corr_file, plot_size,
    divergent_colors, na_color)

  #make cooccur_plot
  plot_file <- file.path(figure_dir, paste("cooccur.", disease_id, ".png",
    sep = ""))
  cooc_plot <- plot_cooccurence(coocur_df, plot_file, plot_size,
    divergent_colors, na_color, q_cut)

  #make disease gene plot
  disease_fig <- file.path(figure_dir, paste("gene_disease.", disease_id,
    ".png", sep = ""))
  disease_plot <- plot_disease(gene_disease_counts, disease_fig, plot_size,
    divergent_colors, na_color)

  #Make combined plot.
  combined_fig <- file.path(figure_dir, paste("combined.", disease_id, ".png",
    sep = ""))
  combine_plots(cooc_plot, disease_plot, combined_fig, coocur_df)

  #write outputs
  cooc_file <- file.path(file_dir, paste("cooccur.", disease_id, ".tsv",
    sep = ""))
  readr::write_tsv(cooccur_summary, cooc_file)
  gene_disease_file <- file.path(file_dir, paste("gene_disease.", disease_id,
    ".tsv", sep = ""))
  readr::write_tsv(gene_disease_counts, gene_disease_file)
  }

print("done")
