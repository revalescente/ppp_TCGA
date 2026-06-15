suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(spatstat))
source("~/ppp_TCGA/src/utils/ppp_fun.R")

option_list <- list(
  make_option(c("--sample_dir"), type = "character", help = "Sample directory"),
  make_option(c("--save_dir"), type = "character", help = "Output directory"),
  make_option(c("--AppMag"), type = "numeric", help = "Magnification factor"),
  make_option(c("--nthreads"), type = "numeric", help = "Thread number for fastKmeans")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Now you can access them like this:
# save_dir <- opt$save_dir
# mask_dir <- opt$mask_dir

# save_path
file_name <- sub("\\.h5ad\\.gz$", ".rds", basename(opt$sample_dir))
save_path <- file.path(opt$save_dir, file_name)
# print(save_path)

# mask loading
mask_path <- sub("/h5ad/", "/mask/", opt$sample_dir)
mask_path <- sub("\\.h5ad\\.gz$", ".png", mask_path)
img <- imager::load.image(mask_path)

# AppMag value is 40 for every one (apparently)

# Convert h5ad to ppp object
pp_object <- ann2ppp(opt$sample_dir, mask = img, magnif = opt$AppMag) # with the image to convert into a mask
# pp_object <- ann2ppp(sample_path, mask = NULL) # without the image but with estimated mask (simple convex hullo)

# creation of vector of distances 
r_values <- seq(0, 150, length.out = 70)
# basic function overflow
#Lcross_neoStroma <- spatstat.explore::Lcross(pp_object, i = "neoplastic", j = "stromal", correction = "border", r = r_values)

# Downsampling random
# Lcross_neoStroma <- Lcross_ds(pp_object, mark_i = "neoplastic", mark_j = "stromal", 
#                               downsample = TRUE, n_down = 40000,
#                               correction = "border", r = r_values)

# Sketch downsampling
# Lcross_neoStroma <- Lcross_sketch(
#   pp_object, 
#   mark_i = "neoplastic", 
#   mark_j = "stromal", 
#   # mark_column = "merged_type", 
#   downsample = TRUE, 
#   n_threads = opt$nthreads,
#   n_down = 10000,
#   correction = "border", 
#   r = r_values
# )

# Lcross with KDE downsampling
set.seed(10)
ppp_sub <- downsample_KDE(pp_object, c("neoplastic", "stromal"), 40000)

Lcross_neoStroma <-spatstat.explore::Lcross(ppp_sub, 
                                            i = "neoplastic", 
                                            j = "stromal", 
                                            r = r_values,
                                            correction = "border")


# saving 
saveRDS(Lcross_neoStroma, file = save_path)


