# 1. Define your parameters
dir_path <- "/projects/shared/TCGA_data/h5ad/"  # The directory to search
n_samples <- 5                                  # The number of samples to select

# 2. List all files in the directory
# full.names = TRUE gives the full path relative to your dir_path
all_files <- list.files(path = dir_path, full.names = TRUE)

# Optional: If you only want specific file types (like .json), use the pattern argument:
# all_files <- list.files(path = dir_path, pattern = "\\.json(\\.gz)?$", full.names = TRUE)

# 3. Safety check: ensure we don't try to sample more files than exist
if (length(all_files) < n_samples) {
  warning("Requested more samples than available files. Selecting all available files.")
  n_samples <- length(all_files)
}

# 4. Randomly select the desired number of files
sampled_files <- sample(all_files, size = n_samples, replace = FALSE)

# 5. Convert to guaranteed absolute paths 
# (Useful in case your initial dir_path was relative, like "./my_folder")
absolute_paths <- normalizePath(sampled_files, mustWork = TRUE)

# Print the result
print(absolute_paths)
write.table(absolute_paths, file = "/projects/shared/TCGA_data/test_samples.txt")

# all samples saved
out_txt <- "/projects/shared/TCGA_data/all_samples.txt"
writeLines(all_files, con = out_txt)
