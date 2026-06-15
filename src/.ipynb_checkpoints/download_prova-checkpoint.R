library(imageFeatureTCGA)
library(purrr)
library(curl)

# 1. Fetch your URLs (as you already did)
url_h5ad    <- getCatalog("hovernet", format = "h5ad") |> dplyr::slice(1) |> getFileURLs()

dir_path <- "/projects/shared/TCGA_data"

# 3. Define the download function with the base directory argument
download_format_files <- function(url_list, format_name, base_dir) {
  
  # Create the full path for the subfolder (e.g., "base_dir/h5ad")
  sub_dir <- file.path(base_dir, format_name)
  
  # Create the directory (recursive = TRUE creates the parent dir_path too if needed)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  
  # Download each file into the new subfolder
  walk(url_list, function(url) {
    file_name <- basename(url)
    dest_path <- file.path(sub_dir, file_name)
    
    if (!file.exists(dest_path)) {
      message("Downloading: ", file_name, " into ", sub_dir)
      curl::curl_download(url, destfile = dest_path, quiet = FALSE)
    } else {
      message("Skipping (already exists): ", file_name)
    }
  })
}

# 4. Run the function for each format, passing the dir_path
download_format_files(url_h5ad, "h5ad", dir_path)

