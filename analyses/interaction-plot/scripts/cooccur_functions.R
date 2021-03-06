# cooccur_function.R
# Functions for calculating co-occurence between mutations

#' Calculate Fisher's exact test for row-wise data
#'
#' Data order follows a two by two matrix, filled by rows or columns
#' (as these are equivalent). Not vectorized!
#'
#' @param w row 1 column 1 value
#' @param x row 2 column 1 value
#' @param y row 1 column 2 value
#' @param z row 2 column 2 value
#'
#' @return p-value from a Fisher's exact test
row_fisher <- function(w, x, y, z) {
  # function to calculate fisher test from row elements
  mat <- matrix(c(w, x, y, z), ncol = 2)
  fisher <- fisher.test(mat)
  return(fisher$p.value)
}

filter_mutations <- function(maf_df, include_var_class) {
  #remove unwanted variant types from the maf_df
  maf_df <- maf_df %>%
    dplyr::filter(Variant_Classification %in% include_var_class)
  return(maf_df)
}

#' Calculate co-occurence relationships for a list of genes
#'
#' @param gene_sample_df a data frame with columns `gene`, `sample`, and
#'   `mutations` where each row represents a gene mutated in the named
#'   sample, and `mutations` is the number of mutations in that gene-sample
#'   combination.
#' @param genes a vector of genes for which co-occurence data is to be calculated.
#'   Default is to pick a random subset, though this is almost never what
#'   is wanted! (Note that all-by-all comparisons is likely to be extremely
#'   slow).
#' @param samples Which samples should the comparison include. Defaults to all
#'   samples present in the data frame.
#'
#' @return A data frame summarizing the co-occurence of pairs of genes in the
#'   gene list with columns `gene1`; `gene2`; counts of each mutations in
#'   each category of sharing (`mut11`: both mutated; `mut10`: mutated in
#'   the first but not second gene, etc.); `odds_ratio` for co-occurence,
#'   `cooccur_sign`: 1 if co-occurence greater than by chance, -1 if less
#'   frequent than expected; `p` the fisher's exact test p value; and a
#'   `cooccur_score` calucated as `cooccur_sign * -log10(p)`.
coocurrence <- function(gene_sample_df,
                        genes = sample(unique(gene_sample_df$gene), 25),
                        samples = unique(gene_sample_df$sample)) {
  # gene_list in order of interest
  #
  # get all pairs of genes
  gene_pairs <- t(combn(genes, m = 2))
  colnames(gene_pairs) <- c("gene1", "gene2")

  # get mutation counts for all genes/sample pairs in gene list
  # fills in any missing values with 0
  all_sample_counts <- expand.grid(
    gene = genes,
    sample = samples,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(gene_sample_df, by = c("gene", "sample")) %>%
    tidyr::replace_na(list(mutations = 0))

  gene_pair_counts <- gene_pairs %>%
    tibble::as_tibble() %>%
    dplyr::left_join(all_sample_counts,
      by = c("gene1" = "gene")
    ) %>%
    dplyr::rename(muts1 = mutations) %>%
    dplyr::left_join(all_sample_counts,
      by = c(
        "gene2" = "gene",
        "sample" = "sample"
      )
    ) %>%
    dplyr::rename(muts2 = mutations) %>%
    dplyr::mutate(
      gene1 = factor(gene1, levels = genes),
      gene2 = factor(gene2, levels = genes)
    )

  gene_pair_summary <- gene_pair_counts %>%
    dplyr::group_by(gene1, gene2) %>%
    dplyr::summarize(
      mut11 = sum(muts1 > 0 & muts2 > 0),
      mut10 = sum(muts1 > 0 & muts2 == 0),
      mut01 = sum(muts1 == 0 & muts2 > 0),
      mut00 = sum(muts1 == 0 & muts2 == 0),
      odds_ratio = (mut11 * mut00) / (mut10 * mut01),
      cooccur_sign = ifelse(odds_ratio > 1, 1, -1)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p = row_fisher(mut11, mut10, mut01, mut00)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      q = p.adjust(p, method = "BH"),
      cooccur_score = cooccur_sign * -log10(p)
    )

  return(gene_pair_summary)
}

modify_cooc_sum <- function(cooccur_summary) {
  #reduce cooccur_summary to necessary fields and generate labels
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
  return(cooccur_df)
}

sig_filter <- function(cooccur_df, p_value) {
  #Remove columns from the cooccur_df with a p_value higher than the cut-off
  filtered_df <- cooccur_df %>%
    dplyr::filter(
      p <= p_value
    )
  return(filtered_df)
}
