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
