file_list <- list.files("/projects/shared/TCGA_data/h5ad", pattern = "\\.h5ad.gz$", full.names = TRUE)
(sample_names <- tools::file_path_sans_ext(basename(file_list)))
(patient_ids_from_files <- substr(sample_names, 1, 12))

meta_matched <- final_metadata |> 
  filter(submitter_id %in% patient_ids_from_files) |> 
  mutate(
    time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    event = ifelse(!is.na(days_to_death), 1, 0)
  )

head(meta_matched)
summary(meta_matched$age_at_diagnosis)

32872/365
meta_matched[which(meta_matched$age_at_diagnosis==32872),]
