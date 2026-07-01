library(ExperimentHub)
library(curatedTCGAData)
library(data.table)

?curatedTCGAData

ov_meta_full <- curatedTCGAData(diseaseCode = "OV", assays = "*", 
                           version = "1.1.38", dry.run = FALSE)

meta_full <- curatedTCGAData(diseaseCode = "*", assays = "*", 
                           version = "1.1.38", dry.run = FALSE)

dt <- as.data.table(colData(ov_meta_full))

colnames(dt)
table(dt$patient.gender)


# ----
# Seconda strada TCGAbiolinks
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(dplyr)
library(purrr)

# Your target disease codes
tcga_disease_codes <- c(
  "ACC",  "BLCA", "BRCA", "CESC", "CHOL", "CNTL", "COAD", "DLBC", 
  "ESCA", "GBM",  "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG",  
  "LIHC", "LUAD", "LUSC", "MESO", "OV",   "PAAD", "PCPG", "PRAD", 
  "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", 
  "UCS",  "UVM"
)

# Function to fetch ONLY clinical data for a given disease code
fetch_clinical_metadata <- function(disease) {
  project_id <- paste0("TCGA-", disease)
  message("Fetching clinical data for: ", project_id)
  
  # Fetch from GDC
  tryCatch({
    data <- GDCquery_clinic(project = project_id, type = "clinical")
    
    # Keep track of the cancer type
    if(nrow(data) > 0) {
      data$disease_code <- disease 
    }
    return(data)
  }, error = function(e) {
    message("Error fetching ", project_id, ": ", e$message)
    return(NULL)
  })
}

# 1. Download and combine all clinical data into one large data frame
all_clinical_data <- map_df(tcga_disease_codes, fetch_clinical_metadata)

# 2. Map your columns to the modern GDC standard names
columns_filtered <- c(
  "submitter_id",               # equivalent to "patient_id"
  "disease_code",               # added by us to track cancer type
  "vital_status",               # "vital_status"
  "year_of_birth",              # proxy for "years_to_birth"
  "days_to_death",              # "days_to_death"
  "days_to_last_follow_up",     # "days_to_last_followup"
  "age_at_diagnosis",           # "patient.age_at_initial_pathologic_diagnosis" (Note: GDC returns this in days)
  "gender",                     # "gender"
  "ethnicity",                  # "ethnicity"
  "tissue_or_organ_of_origin",   # "patient.anatomic_neoplasm_subdivision"
  "race",                   # Essential demographic pairing with ethnicity
  "ajcc_pathologic_stage",  # Tumor stage
  "tumor_stage",            # Sometimes stage is recorded here instead depending on the cancer
  "primary_diagnosis",      # Histological subtype
  "prior_malignancy",        # Useful for filtering out confounding previous cancers  
  "tissue_source_site" 
)

# 3. Filter the dataframe to only keep the columns of interest 
# (intersect ensures we don't crash if a column happens to be missing)
existing_cols <- intersect(columns_filtered, colnames(all_clinical_data))
final_metadata <- all_clinical_data %>% 
  select(all_of(existing_cols)) %>%
  # 3. Create a reliable 'hospital_code' column directly from the ID
  # The code is always the 6th and 7th character in "TCGA-XX-YYYY"
  mutate(hospital_code = substr(submitter_id, 6, 7))

# View the result
colnames(final_metadata)
head(final_metadata)
# View the result
head(final_metadata)

