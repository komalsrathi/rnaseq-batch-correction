#R script wrapper for generating interaction plots

#TODO / Questions:
# Add CNV options and parsing?
# Flags for different mutation types (syn, non-coding)

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
    help = "File path to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outdir",
    default = file.path(getwd(), "out"),
    type = "character",
    help="Output directory path."
  ),
  make_option(
    opt_str = "--exclude",
    default = file.path(getwd(), "exclude-genes.txt"),
    type = "character",
    help = "File path with a table of genes to be excluded from the figure.
      A tsv file which must contain a column named 'gene` that contains Hugo Symbols"
  ),
  make_option(
    opt_str = "--samples",
    type = "character",
    help = "File path to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--maf",
    type = "character",
    help = "File path of MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--exclude_genes",
    type = "character",
    default = NA,

  ),
  make_option(
    opt_str = "--scripts",
    default = file.path(getwd(), "scripts"),
    type = "character",
    help="Scripts directory path."
  )
)

#defaults that might be changed later
max_genes <- 50
min_vaf <- 0.05 #minimum vaf for inclusion
min_mutated <- 5 #Number of samples that have to have mutation for it to be included
min_depth <- 0 #minimum depth of coverage

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
meta_file <- file.path(opts$metadata)
sample_file <- file.path(opts$samples)
maf_file <- file.path(opts$maf)
out_dir <- file.path(opts$outdir)
exclude_file <- file.path(opts$exclude)
script_dir <- file.path(opts$scripts)
print(meta_file)
print(out_dir)
print(script_dir)

#load additional functions
#we might only need the one file, so it might be better to take
#the path to that script not the folder
source(file.path(script_dir, "cooccur_functions.R"))

#create output directories
figure_dir <- file.path(out_dir, "figures")
file_dir <- file.path(out_dir, "files")
dir.create(out_dir)
dir.create(figure_dir)
dir.create(file_dir)

#set up diseases
diseases <- hash(
"all"="All",
"Medulloblastoma"="Medulloblastoma",
"LGAT"="Low-grade astrocytic tumor",
"Ependymoma"="Ependymoma",
"HGAT"="High-grade glioma",
"DMG"="Diffuse midline glioma",
"Ganglioglioma"="Ganglioglioma",
"Craniopharyngioma"="Craniopharyngioma"
)

#read inputs
meta_df <- readr::read_tsv(meta_file, col_types = readr::cols(), guess_max = 10000)
sample_df <- readr::read_tsv(sample_file, col_types = readr::cols())
maf_df <- data.table::fread(maf_file, data.table = FALSE)
exclude_df <- data.table::fread(exclude_file, data.table = FALSE)

### Reduce MAF to a smaller set of relevant columns
maf_df <- maf_df %>%
  dplyr::select(
    Hugo_Symbol,
    Entrez_Gene_Id,
    Chromosome,
    Start_Position,
    End_Position,
    Strand,
    Variant_Classification,
    Variant_Type,
    Reference_Allele,
    Tumor_Seq_Allele1,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    t_depth,
    t_ref_count,
    t_alt_count,
    Consequence
  )

#get ncbi file (or take as option?)
#get gene list
genes <- unique(maf_df$Hugo_Symbol)
genes <- genes[!genes %in% exclude_df$gene] #remove genes from exclude list

#organize variant types and decide which to analyze
# generate consequence lists and filter
intergenic <- c("IGR")
nontranscribed <- c(
  "3'Flank",
  "5'Flank",
  "Targeted_Region"
)
noncoding <- c(
  "RNA",
  "Intron",
  "3'UTR",
  "5'UTR",
  "Splice_Region",
  "lincRNA"
)
# Variant Classification with Low/Modifier variant consequences
#  from maftools http://asia.ensembl.org/Help/Glossary?id=535
synonymous <- c(
  "Silent",
  "Start_Codon_Ins",
  "Start_Codon_SNP",
  "Stop_Codon_Del",
  "De_novo_Start_InFrame",
  "De_novo_Start_OutOfFrame"
)
# Variant Classification with High/Moderate variant consequences from maftools
nonsynonymous <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)
include <- nonsynonymous # always want nonsyn
#if (opts$include_syn) {
#  include <- c(include, synonymous)
#}
#if (opts$include_noncoding) {
#  include <- c(include, noncoding)
#}
#if (opts$include_nontranscribed) {
#  include <- c(include, nontranscribed)
#}

#run analysis for each disease
for (diseaseID in keys(diseases))
  {
  print(diseaseID)
  full_diseaseID <- diseases[[diseaseID]]
  print(full_diseaseID)
  #filter for samples with the disease
  disease_df <- meta_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% sample_df$Kids_First_Biospecimen_ID) %>%
    dplyr::filter(tolower(full_diseaseID) == "all" |
      tolower(integrated_diagnosis) == tolower(full_diseaseID)) %>%
    dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

  #filter metadata to chosen samples
  samples <- disease_df$Kids_First_Biospecimen_ID
  sample_meta <- meta_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% samples)

  #reduce maf to chosen samples & calculate VAF
  sample_maf <- maf_df %>%
  dplyr::filter(Tumor_Sample_Barcode %in% samples) %>%
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count))
  #filter maf file
  maf_filtered <- sample_maf %>%
  filter_mutations(
    min_vaf = min_vaf,
    min_depth = min_depth,
    include_var_class = include
  )

  # count mutations by gene/sample pair
  gene_sample_counts <- maf_filtered %>%
  dplyr::filter(Entrez_Gene_Id > 0, # remove unknowns
    Hugo_Symbol %in% genes) %>% # include only desired genes
  dplyr::group_by(gene = Hugo_Symbol, sample = Tumor_Sample_Barcode) %>%
  dplyr::tally(name = "mutations") %>%
  dplyr::ungroup()
  # count # of samples mutated by gene
  gene_counts <- gene_sample_counts %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(
    mutant_samples = dplyr::n(),
    total_muts = sum(mutations),
    mean_muts_per_sample = mean(mutations)
    ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    desc(mutant_samples),
    desc(mean_muts_per_sample)
    ) %>%
  dplyr::filter(mutant_samples >= min_mutated |
  dplyr::row_number() <= 2) # keep at least 2 genes

  #get most often mutated genes
  top_count_genes <- head(gene_counts, max_genes)$gene
  cooccur_summary <- coocurrence(gene_sample_counts, top_count_genes)

  #count mutated samples by disease type
  gene_disease_counts <- gene_sample_counts %>%
  dplyr::filter(gene %in% top_count_genes) %>%
  dplyr::left_join(sample_meta,
    by = c("sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::group_by(gene, disease = integrated_diagnosis) %>%
  dplyr::summarize(mutant_samples = dplyr::n(),
    total_muts = sum(mutations),
    mean_muts_per_sample = mean(mutations)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    desc(mutant_samples),
    desc(mean_muts_per_sample)
  )

  #write outputs
  disease_file <- file.path(file_dir, paste(diseaseID, ".tsv", sep=''))
  readr::write_tsv(disease_df, disease_file)
  cooc_file <- file.path(file_dir, paste("cooccur.", diseaseID, ".tsv", sep=''))
  readr::write_tsv(cooccur_summary, cooc_file)
  gene_disease_file <- file.path(file_dir, paste("gene_disease.", diseaseID, ".tsv", sep=''))
  readr::write_tsv(gene_disease_counts, gene_disease_file)
  }

print("done")
