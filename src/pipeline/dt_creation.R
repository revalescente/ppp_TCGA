suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

# Directories
input_dir <- "/projects/shared/TCGA_data/Lfun_KDE"
# Changed extension to .csv
output_file <- "/projects/shared/TCGA_data/combined_Lcross_KDE.csv" 

# 1. Get files and names
files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
sample_names <- sub("\\.rds$", "", basename(files))

# 2. Get the number of cores from SLURM environment (default to 4 if not found)
slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
n_cores <- if (slurm_cores != "") as.numeric(slurm_cores) else 4

cat(sprintf("Found %d files.\n", length(files)))
cat(sprintf("Starting parallel extraction using %d CPU cores...\n", n_cores))

# 3. Read and extract in parallel
res_list <- mclapply(files, function(f) {
  tryCatch({
    obj <- readRDS(f)
    return(obj$border)
  }, error = function(e) {
    # If a file is corrupted, return NAs to keep the row dimensions matching
    return(rep(NA_real_, 70)) 
  })
}, mc.cores = n_cores)

# 4. Name list and convert to datatable
names(res_list) <- sample_names
dt <- as.data.table(res_list)

# 5. Extract r values from the first valid file and prepend
first_obj <- readRDS(files[1])
dt <- data.table(r_value = first_obj$r, dt)

# 6. Save final output using fwrite (lightning fast CSV writer)
fwrite(dt, file = output_file)

cat("\nDimensions of final table:", nrow(dt), "rows by", ncol(dt), "columns.\n")
cat("Successfully saved to:", output_file, "\n")