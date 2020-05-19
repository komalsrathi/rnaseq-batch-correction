#R package to organize samples and genes by disease type
#and calculate coocurrence score

library(dplyr)

samples_with_disease <- function(meta_df, sample_df, disease_id) {
  #returns a list of samples with given disease from the meta_df
  disease_df <- meta_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% sample_df$Kids_First_Biospecimen_ID) %>%
    dplyr::filter(tolower(disease_id) == "all" |
      tolower(short_histology) == tolower(disease_id)) %>%
    dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

  samples <- disease_df$Kids_First_Biospecimen_ID
  return(samples)
}

gene_counts <- function (maf_df, genes) {
  #count the number of samples with a mutation in each gene
  gene_sample_counts <- maf_df %>%
  dplyr::filter(Entrez_Gene_Id > 0, # remove unknowns
    Hugo_Symbol %in% genes) %>% # include only desired genes
    dplyr::group_by(gene = Hugo_Symbol, sample = Tumor_Sample_Barcode) %>%
    dplyr::tally(name = "mutations") %>%
    dplyr::ungroup()
  return(gene_sample_counts)
}

get_top_genes <- function (gene_sample_counts, max_genes, min_mutated) {
  # count samples with mutations in each gene
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

  return(top_count_genes)
}
