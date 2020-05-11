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

get_gene_list <- function(maf_df, exclude_file, disease_id) {
  #given the maf file and disease id, return gene list
  #only exclude commonly mutated genes when the disease group is all
  #get ncbi file (or take as option?) for exclude genes
  #get initial gene list
  genes <- unique(maf_df$Hugo_Symbol)
  if (tolower(disease_id) == "all") {
    exclude_df <- data.table::fread(exclude_file, data.table = FALSE)
    genes <- genes[!genes %in% exclude_df$gene] #remove genes from exclude list
  }
  return(genes)
}

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
      A tsv file which must contain a column named 'gene` that contains Hugo Symbols"
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
  )
)

#defaults that might be changed later
max_genes <- 50
min_vaf <- 0.05 #minimum vaf for inclusion
min_mutated <- 5 #Number of samples that have to have mutation for it to be included
min_depth <- 0 #minimum depth of coverage
plot_size <- 50

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
meta_file <- file.path(opts$metadata)
sample_file <- file.path(opts$samples)
maf_file <- file.path(opts$maf)
out_dir <- file.path(opts$outdir)
exclude_file <- file.path(opts$exclude)
script_dir <- file.path(opts$scripts)
pallette_file <- file.path(opts$palette)
print(meta_file)
print(out_dir)
print(script_dir)

#load additional functions
#we might only need the one file, so it might be better to take
#the path to that script not the folder
source(file.path(script_dir, "cooccur_functions.R"))

#set up color palette
divergent_palette <- readr::read_tsv(pallette_file, col_types = readr::cols())
divergent_colors <- divergent_palette %>%
  dplyr::filter(color_names != "na_color") %>%
  dplyr::pull(hex_codes)
na_color <- divergent_palette %>%
  dplyr::filter(color_names == "na_color") %>%
  dplyr::pull(hex_codes)

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

#read inputs
meta_df <- readr::read_tsv(meta_file, col_types = readr::cols(), guess_max = 10000)
sample_df <- readr::read_tsv(sample_file, col_types = readr::cols())
maf_df <- data.table::fread(maf_file, data.table = FALSE)

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
for (disease_id in keys(diseases)) {
  print(disease_id)
  full_disease_id <- diseases[[disease_id]]
  print(full_disease_id)

  genes <- get_gene_list(maf_df, exclude_file, disease_id)

  #filter for samples with the disease
  disease_df <- meta_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% sample_df$Kids_First_Biospecimen_ID) %>%
    dplyr::filter(tolower(full_disease_id) == "all" |
      tolower(integrated_diagnosis) == tolower(full_disease_id)) %>%
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

  #reduce cooccur summary to needed fields
  cooccur_df <- cooccur_summary %>%
    dplyr::mutate(
    mut1 = mut11 + mut10,
    mut2 = mut11 + mut01,
    label1 = paste0(gene1, " (", mut1, ")"),
    label2 = paste0(gene2, " (", mut2, ")")
  )

  labels <- unique(c(cooccur_df$label1, cooccur_df$label2))

  #check the order of the labels to be decreasing by mut count
  label_counts <- as.numeric(stringr::str_extract(labels, "\\b\\d+\\b"))
  labels <- labels[order(label_counts, decreasing = TRUE)]
  #order genes the same way, in case we want to use those
  genes <- stringr::str_extract(labels, "^.+?\\b")
  genes <- genes[order(label_counts, decreasing = TRUE)]

  cooccur_df <- cooccur_df %>%
    dplyr::mutate(
    gene1 = factor(gene1, levels = genes),
    gene2 = factor(gene2, levels = genes),
    label1 = factor(label1, levels = labels),
    label2 = factor(label2, levels = labels)
  )

  #create scales for consistent sizing
  #The scales will need to have plot_size elements,
  #so after getting the unique list, we concatenate on extra elements.
  #for convenience, these are just numbers 1:n
  #where n is the number of extra labels needed for the scale
  xscale <- cooccur_df$label1 %>%
    as.character() %>%
    unique() %>%
    c(1:(plot_size - length(.)))
  yscale <- cooccur_df$label2 %>%
    as.character() %>%
    unique() %>%
    #the concatenated labels need to be at the front of the Y scale,
    #since this will be at the bottom in the plot.
    c(1:(plot_size - length(.)), .)

    #make cooccur plot
    cooccur_plot <- ggplot(
      cooccur_df,
      aes(x = label1, y = label2, fill = cooccur_score)
    ) +
      geom_tile(width = 0.7, height = 0.7) +
      scale_x_discrete(
        position = "top",
        limits = xscale,
        breaks = unique(cooccur_df$label1)
      ) + # blank unused sections.
      scale_y_discrete(
        limits = yscale,
        breaks = unique(cooccur_df$label2)
      ) +
      scale_fill_gradientn(
        colors = divergent_colors,
        na.value = na_color,
        limits = c(-10, 10),
        oob = scales::squish,
      ) +
      labs(
        x = "",
        y = "",
        fill = "Co-occurence\nscore"
      ) +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        axis.text.x = element_text(
          angle = -90,
          hjust = 1,
          size = 6
        ),
        axis.text.y = element_text(size = 6),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.key.size = unit(2, "char")
      )
  plot_file <- file.path(figure_dir, paste("cooccur.", disease_id, ".png", sep = ""))
  print(plot_file)
  ggsave(cooccur_plot, filename = plot_file)

  #Make gene by disease chart.
  #this part needs to be adapted and changed.
  #the original code is calling the gene_disease_counts disease_df
  #this new code is using a different disease_df
  #disease_df_save <- disease_df #save full disease_df for printing later.
  disease_gene_df <- gene_disease_counts %>%
    dplyr::mutate(gene = factor(gene, levels = genes))
  display_diseases <- disease_gene_df %>%
    dplyr::group_by(disease) %>%
    dplyr::tally(wt = mutant_samples) %>%
    dplyr::arrange(desc(n)) %>%
    head(7) %>% # seven so we end up with 8 total for color reasons
    dplyr::pull(disease)
  disease_gene_df <- disease_gene_df %>%
    dplyr::mutate(disease_factor =
      forcats::fct_other(disease, keep = display_diseases) %>%
      forcats::fct_relevel(display_diseases)
    )
  xscale2 <- levels(disease_gene_df$gene) %>%
      c(rep("", plot_size - length(.)))
  disease_plot <- ggplot(
    disease_gene_df,
    aes(x = gene, y = mutant_samples, fill = disease_factor)) +
    geom_col(width = 0.7) +
    labs(
      x = "",
      y = "Samples with mutations",
      fill = "Diagnosis"
      ) +
    colorblindr::scale_fill_OkabeIto() +
    scale_x_discrete(
      limits = xscale2,
      breaks = disease_gene_df$gene
    ) +
    scale_y_continuous(expand = c(0, 0.5, 0.1, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.key.size = unit(1, "char"))
  #save disease plot
  disease_fig <- file.path(figure_dir, paste("gene_disease.", disease_id, ".png", sep = ""))
  print(disease_fig)
  ggsave(disease_plot, filename = disease_fig)

  #Make combined plot.
  #labels for y axis will be gene names, with extra spaces (at bottom) blank
  ylabels  <- cooccur_df$gene2 %>%
    as.character() %>%
    unique() %>%
    c(rep("", plot_size - length(.)), .)
  cooccur_plot2 <- cooccur_plot +
    scale_x_discrete(
      limits = xscale,
      breaks = c()
    ) +
    scale_y_discrete(
      limits = yscale,
      labels = ylabels
    ) +
    theme(
      plot.margin = unit(c(-3.5, 0, 0, 0), "char")# negative top margin to move plots together
    )
  #Move labels and themes for disease plot
  disease_plot2 <- disease_plot +
    theme(
      axis.text.x = element_text(
        angle = -90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.title.y = element_text(
        vjust = -10# keep the label close when combined
      )
    )

  #Combine plots with <patchwork>
  #Layout of the two plots will be one over the other (1 column),
  #with the upper plot 3/4 the height of the lower plot
  combined_plot <- disease_plot2 + cooccur_plot2 +
    plot_layout(ncol = 1, heights = c(3, 4)) +
    plot_annotation(tag_levels = "A") &
    theme(#add uniform labels
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )
  #save combined plot.
  combined_fig <- file.path(figure_dir, paste("combined.", disease_id, ".png", sep = ""))
  print(combined_fig)
  ggsave(combined_plot, filename = combined_fig, width = 8, height = 14)

  #write outputs
  disease_file <- file.path(file_dir, paste(disease_id, ".tsv", sep = ""))
  readr::write_tsv(disease_df, disease_file)
  cooc_file <- file.path(file_dir, paste("cooccur.", disease_id, ".tsv", sep = ""))
  readr::write_tsv(cooccur_summary, cooc_file)
  gene_disease_file <- file.path(file_dir, paste("gene_disease.", disease_id, ".tsv", sep = ""))
  readr::write_tsv(gene_disease_counts, gene_disease_file)
  }

print("done")
