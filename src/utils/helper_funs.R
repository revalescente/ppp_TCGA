library(purrr)
library(data.table)

extract_tcga_clinical <- function(disease_codes, columns_to_select, version = "2.1.1") {
  
  # Usiamo map (restituisce una lista di data.table) invece di map_dfr
  list_of_dts <- disease_codes |> 
    purrr::map(\(code) {
      message("---> Elaborazione del cancro: ", code)
      
      tryCatch({
        # 1. Scarica il MultiAssayExperiment minimo
        mae <- curatedTCGAData(diseaseCode = code, 
                               assays = "RNASeqGene", 
                               version = version, 
                               dry.run = FALSE)
        
        # 2. Converti in data.table. 
        # NOTA: data.table rimuove i rownames, quindi usiamo keep.rownames 
        # per salvare l'ID del campione in una colonna reale chiamata "sample_id"
        dt <- as.data.table(as.data.frame(colData(mae)), keep.rownames = "sample_id")
        
        # 3. Seleziona solo le colonne richieste che esistono in questo tumore
        existing_cols <- intersect(columns_to_select, names(dt))
        dt <- dt[, c("sample_id", existing_cols), with = FALSE]
        
        # 4. Aggiungi il codice del tumore alla velocità della luce (in-place)
        dt[, disease_code := code]
        
        return(dt)
        
      }, error = function(e) {
        message("⚠️ Errore con il codice ", code, ": ", e$message)
        return(NULL) # Verrà ignorato da rbindlist
      })
    })
  
  # Il superpotere di data.table: unisce la lista di tabelle istantaneamente
  # e inserisce gli NA dove un tumore non ha una determinata colonna.
  final_dataset <- data.table::rbindlist(list_of_dts, fill = TRUE)
  
  return(final_dataset)
}
