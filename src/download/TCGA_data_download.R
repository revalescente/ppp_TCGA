library(purrr)
library(curl)

# 1. Define the main directory where everything should be saved
dir_path <- "/projects/shared/TCGA_data" # <-- CHANGE THIS to your desired path

# 2. Fetch your URLs
url_h5ad    <- readRDS('~/ppp_TCGA/data/catalogs/h5ad_TCGAcatalog.rds')
url_geojson <- readRDS('~/ppp_TCGA/data/catalogs/geojson_TCGAcatalog.rds') 
url_thumb   <- readRDS('~/ppp_TCGA/data/catalogs/thumb_TCGAcatalog.rds') 
url_json    <- readRDS('~/ppp_TCGA/data/catalogs/json_TCGAcatalog.rds') 


# 3. Define the robust download function with error handling
download_format_files <- function(url_list, format_name, base_dir) {
  
  sub_dir <- file.path(base_dir, format_name)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  walk(url_list, function(url) {
    file_name <- basename(url)
    dest_path <- file.path(sub_dir, file_name)
    
    if (!file.exists(dest_path)) {
      message("Downloading: ", file_name, " into ", sub_dir)
      
      # tryCatch prevents the script from crashing on a network error
      tryCatch({
        curl::curl_download(url, destfile = dest_path, quiet = FALSE)
      }, error = function(e) {
        # If it fails, print the error
        message("\n[ERROR] Failed to download ", file_name, ": ", e$message)
        
        # CRITICAL: Delete the partial file so it isn't skipped on the next run!
        if (file.exists(dest_path)) {
          file.remove(dest_path)
          message("Cleaned up partial file for ", file_name)
        }
      })
      
    } else {
      message("Skipping (already exists): ", file_name)
    }
  })
}

# 4. Run the function for each format, passing the dir_path
download_format_files(url_h5ad, "h5ad", dir_path)
download_format_files(url_geojson, "geojson", dir_path)
download_format_files(url_thumb, "thumb", dir_path)
download_format_files(url_json, "json", dir_path)



