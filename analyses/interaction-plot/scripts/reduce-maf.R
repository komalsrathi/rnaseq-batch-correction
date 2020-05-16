#R package to reduce maf data frame to needed fields
#will probably be quickly deprecated

library(dplyr)

reduce_maf <- function (maf_df) {
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
  return(maf_df)
}
