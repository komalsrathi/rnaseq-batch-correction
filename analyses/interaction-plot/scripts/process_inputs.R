get_gene_list <- function(maf_df, exclude_file, disease_id) {
  #given the maf data.frame and disease id, return gene list
  #only exclude commonly mutated genes when the disease group is all
  #get initial gene list
  genes <- unique(maf_df$Hugo_Symbol)
  if (tolower(disease_id) == "all") {
    exclude_df <- data.table::fread(exclude_file, data.table = FALSE)
    genes <- genes[!genes %in% exclude_df$gene] #remove genes from exclude list
  }
  return(genes)
}
