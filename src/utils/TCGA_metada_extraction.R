# lista dei tipi di cancro
# 1     ACC                                          Adrenocortical Carcinoma
# 2    BLCA                                      Bladder Urothelial Carcinoma
# 3    BRCA                                         Breast Invasive Carcinoma
# 4    CESC  Cervical Squamous Cell Carcinoma And Endocervical Adenocarcinoma
# 5    CHOL                                                Cholangiocarcinoma
# 6    CNTL                                                          Controls
# 7    COAD                                              Colon Adenocarcinoma
# 8    DLBC                   Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
# 9    ESCA                                              Esophageal Carcinoma
# 10    GBM                                           Glioblastoma Multiforme
# 11   HNSC                             Head And Neck Squamous Cell Carcinoma
# 12   KICH                                                Kidney Chromophobe
# 13   KIRC                                 Kidney Renal Clear Cell Carcinoma
# 14   KIRP                             Kidney Renal Papillary Cell Carcinoma
# 15   LAML                                            Acute Myeloid Leukemia
# 16    LGG                                          Brain Lower Grade Glioma
# 17   LIHC                                    Liver Hepatocellular Carcinoma
# 18   LUAD                                               Lung Adenocarcinoma
# 19   LUSC                                      Lung Squamous Cell Carcinoma
# 20   MESO                                                      Mesothelioma
# 21     OV                                 Ovarian Serous Cystadenocarcinoma
# 22   PAAD                                         Pancreatic Adenocarcinoma
# 23   PCPG                                Pheochromocytoma And Paraganglioma
# 24   PRAD                                           Prostate Adenocarcinoma
# 25   READ                                             Rectum Adenocarcinoma
# 26   SARC                                                           Sarcoma
# 27   SKCM                                           Skin Cutaneous Melanoma
# 28   STAD                                            Stomach Adenocarcinoma
# 29   TGCT                                       Testicular Germ Cell Tumors
# 30   THCA                                                 Thyroid Carcinoma
# 31   THYM                                                           Thymoma
# 32   UCEC                              Uterine Corpus Endometrial Carcinoma
# 33    UCS                                            Uterine Carcinosarcoma
# 34    UVM                                                    Uveal Melanoma
library(ExperimentHub)
library(curatedTCGAData)
ov_meta <- curatedTCGAData(diseaseCode = "OV", assays = "RNASeqGene", version = "2.1.1", dry.run = FALSE)

# prendo questo sottoinsieme di colonne per adesso, per analisi preliminari
columns_to_select <- c(
  "patient_id",                             # patient id
  "vital_status",                          # Survival status
  "years_to_birth",                          # Birth date
  "days_to_death",                         # Time variable for deceased patients
  "days_to_last_followup",                 # Time variable for censored patients
  "patient.age_at_initial_pathologic_diagnosis", # Age 
  "gender",                                # Gender/Sex
  "ethnicity",                             # Ethnicity
  "patient.anatomic_neoplasm_subdivision" # Anatomic site 
  # TBD adding other vars of interest
  
)
colnames(ov_meta)[colnames(ov_meta) %in% columns_to_select]
meta_sub <- ov_meta[, columns_to_select, drop = FALSE]
meta_sub |> head()


library(curatedTCGAData)
library(MultiAssayExperiment)

extract_tcga_clinical <- function(disease_codes, 
                                  columns_to_select, 
                                  version = "2.1.1") {
  
  # Lista temporanea per salvare i data.frame di ogni tumore
  combined_list <- list()
  
  for (code in disease_codes) {
    message("---> Elaborazione del cancro: ", code)
    
    # Gestione degli errori per evitare che il ciclo si blocchi se un download fallisce
    df_clean <- tryCatch({
      
      # 1. Scarica l'oggetto MultiAssayExperiment minimo
      mae <- curatedTCGAData(diseaseCode = code, 
                             assays = "RNASeqGene", 
                             version = version, 
                             dry.run = FALSE)
      
      # 2. Estrai i colData e convertili in un data.frame classico di R
      clinical_df <- as.data.frame(colData(mae))
      
      # 3. Creiamo un data.frame vuoto con le colonne desiderate (salva-errore)
      # Se una colonna non esiste nel tumore X, verrà riempita di NA
      sub_df <- data.frame(matrix(NA, nrow = nrow(clinical_df), ncol = length(columns_to_select)))
      colnames(sub_df) <- columns_to_select
      rownames(sub_df) <- rownames(clinical_df)
      
      # Individua quali colonne richieste sono effettivamente presenti
      existing_cols <- intersect(columns_to_select, colnames(clinical_df))
      
      # Copia solo le colonne esistenti
      if (length(existing_cols) > 0) {
        sub_df[, existing_cols] <- clinical_df[, existing_cols, drop = FALSE]
      }
      
      # 4. Aggiungi il codice del tumore come colonna di tracciamento
      sub_df$disease_code <- code
      
      sub_df
      
    }, error = function(e) {
      message("⚠️ Errore con il codice ", code, ": ", e$message)
      return(NULL)
    })
    
    # Se il download e l'estrazione sono andati a buon fine, salva nella lista
    if (!is.null(df_clean)) {
      combined_list[[code]] <- df_clean
    }
  }
  
  # Unisci tutti i data.frame della lista in un unico grande dataset finale
  if (length(combined_list) > 0) {
    final_dataset <- do.call(rbind, combined_list)
    return(final_dataset)
  } else {
    stop("Nessun dato estratto. Controlla la connessione o i codici inseriti.")
  }
}

# Le tue colonne di interesse
my_columns <- c(
  "patient_id", "vital_status", "years_to_birth",
  "days_to_death", "days_to_last_followup",
  "patient.age_at_initial_pathologic_diagnosis", 
  "gender", "ethnicity", "patient.anatomic_neoplasm_subdivision"
)

# Estrazione singola
ov_clinical_only <- extract_tcga_clinical(disease_codes = "OV", columns_to_select = my_columns)

head(ov_clinical_only)