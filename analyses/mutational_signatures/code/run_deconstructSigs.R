#' Author: C. Savonen for ALSF CCDL and Krutika Gaonkar for D3b
#'
#'
#' Reads in maf file and runs deconstructSigs to get weights per signature  along with muts per signature
#' @param maf maf format mutation file
#' @param genomeversion should be hg19 or hg38[Default]
#' @param normalization deconstructSigs provides the following options : "default"[Default], "exome","genome","exome2genome"
#' @param signatures if "cosmic" [Default] use deconstructSigs provided signatures.cosmic if "nature2013" use deconstructSigs provided signatures.nature2013
#' @param wgs_genome_size: size of the WGS genome in bp
#' @param wxs_genome_size: size of the WXS exome in bp
#' @param metadata_df a data.frame with grouping column and `experimental strategy`  information columns
#' @param grouping_by in the metadata the column name of the grouping variable `short_histology`[Default] 
#' @return A data.frame that is samples x signatures and has the number of mutations per mb for each signature x sample combo



run_deconstructSigs<-function(maf=maf,genomeversion="hg38",normalization="default",signatures="cosmic",wgs_bed=NULL,wxs_bed=NULL,metadata_df=NULL,grouping_by="short_histology"){
  if (genomeversion =="hg19"){
    # Convert maf to deconstructSigs input ;defaults to hg19 Bsgenome
    sigs_input <- mut.to.sigs.input(mut.ref = maf, 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chromosome", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Allele")
    
  }
  if (genomeversion =="hg38"){
    # need BSgenome.Hsapiens.UCSC.hg38 because hg19 is run by default  
    if(!("BSgenome.Hsapiens.UCSC.hg38" %in% installed.packages())){
      BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    }
    
    # Convert maf to deconstructSigs input
    sigs_input <- mut.to.sigs.input(mut.ref = maf, 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chromosome", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Allele",
                                    bsg = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  }
  
  # Set up BED region files for TMB calculations
  if(!is.null(wgs_bed) & !is.null(wxs_bed)){
    wgs_bed <- readr::read_tsv(wgs_bed,col_names = FALSE)
    wxs_bed <- readr::read_tsv(wxs_bed,col_names = FALSE)
    # Calculate size of genome surveyed
    # These files are BED files where the third column is the End position and
    # the second column is the Start position.
    # So End - Start gives the size of each range. Sum the gives the total size in bp.
    wgs_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
    wxs_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])
  } else if(!is.null(wgs_bed) & is.null(wxs_bed)){
    # Calculate size of genome surveyed
    wgs_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
  } else if(is.null(wgs_bed) & !is.null(wxs_bed)){
    # Calculate size of genome surveyed
    wxs_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])
  } else{
    warning("wgs nor wxs bed provided")
  }
  
  # get samples that were retained as deconstructSigs input in sigs_input
  tumor_sample_ids <- maf %>%
    dplyr::filter(Tumor_Sample_Barcode %in% rownames(sigs_input)) %>%
    dplyr::distinct(Tumor_Sample_Barcode) %>%
    dplyr::pull(Tumor_Sample_Barcode)
  
  # Count the total number of signature mutations for each sample
  total_muts <- apply(sigs_input, 1, sum)
  
  # check if signature matrix is available
  if(signatures=="cosmic"){
    signatures.ref<-signatures.cosmic
  } else if (signatures=="nature2013"){
    signatures.ref<-signatures.nature2013
  }  else{
    stop("signatures can be cosmic or nature2013")
  }
  
  if (any(normalization %in% c("default","exome","genome","exome2genome"))){
  sample_sigs <- lapply(tumor_sample_ids, function(sample_id) {
    # Determine the signatures contributing to the sample
    whichSignatures(
      tumor.ref = sigs_input,
      signatures.ref = signatures.ref,
      sample.id = sample_id,
      contexts.needed = TRUE,
      tri.counts.method = normalization,
    )
  })
  # Bring along the names
  names(sample_sigs) <- tumor_sample_ids
  }else{
    stop("Provide normalization methods")
  }
  
  # Pull out the signature weights and make into matrix
  sig_num_df <- do.call(
    "rbind.data.frame",
    lapply(sample_sigs, function(sample_data) sample_data$weights)
  ) %>%
    tibble::rownames_to_column("Tumor_Sample_Barcode") %>%
    
    # Calculate the number of mutations contributing to each signature
    # Here the weight is multiplied by the total number of signature mutations. 
    dplyr::mutate_at(dplyr::vars(-Tumor_Sample_Barcode), ~ . * total_muts) %>%
    
    # Join the grouping variable and experimental stategy information
    dplyr::left_join(metadata_df %>% 
                       dplyr::select(!!(as.name(grouping_by)), experimental_strategy,Tumor_Sample_Barcode) %>%
                       dplyr::distinct(Tumor_Sample_Barcode,.keep_all = TRUE),
                     by = "Tumor_Sample_Barcode"
    ) %>%
    
    # Get rid of Panel samples
    dplyr::filter(experimental_strategy != "Panel") %>%
    
    # Reformat for plotting
    reshape2::melt(value.name = "num_mutations") %>%
    
    # Add genome size and calculate the mutation per this column
    dplyr::mutate(
      genome_size = dplyr::recode(experimental_strategy,
                                  "WGS" = wgs_size,
                                  "WXS" = wxs_size
      ),
      mut_per_mb = num_mutations / (genome_size / 10^6)
    ) %>% 
    # Rename variable as signature
    dplyr::rename("signature" = variable) %>% 
    
    # Make this a factor but make sure the levels stay in order
    dplyr::mutate(signature = factor(signature, levels = unique(signature)))
  
  
  return(sig_num_df)  
}
