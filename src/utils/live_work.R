library(Morpho)
library(spatstat)
library(ggplot2)
library(patchwork)
library(parallel)
library(purrr)
source("~/ppp_TCGA/src/utils/ppp_fun.R")

# campione grande
sample_path <- "/projects/shared/TCGA_data/h5ad/TCGA-4X-A9FB-01Z-00-DX1.211CC9AA-F721-4D16-8663-68A393223F80.h5ad.gz"


# campione problematico
#sample_path <- "/projects/shared/TCGA_data/h5ad/TCGA-06-0213-01Z-00-DX1.16219fe9-0be9-46b2-9190-02b981c4bc39.h5ad.gz"
(file_name <- sub("\\.h5ad\\.gz$", "", basename(sample_path)))


# mask loading
mask_path <- sub("/h5ad/", "/mask/", sample_path)
mask_path <- sub("\\.h5ad\\.gz$", ".png", mask_path)
img <- imager::load.image(mask_path)

# mpp values

# Convert h5ad to ppp object
pp_object2 <- ann2ppp(sample_path, mask = img, magnif = 40)

# sketching 
# extract coord of the specified cell type
ind_i <- which(marks(pp_object) == "neoplastic")
ind_j <- which(marks(pp_object) == "stromal")
coords <- cbind(x = pp_object$x[ind_i], y = pp_object$y[ind_i])
# fastKeamns from Morpho
km <- Morpho::fastKmeans(coords, k = 10000, iter.max = 15, threads = 8)

# Fast Vectorized Sampling
# Randomly shuffle all the indices
rand_order <- sample(seq_along(ind_i))
# Match the clusters to the shuffled order
shuffled_clusters <- km$class[rand_order]
# Pick the first occurrence
first_occurrences <- !duplicated(shuffled_clusters)
# Extract the selected relative indices
sketched_rel_idx <- rand_order[first_occurrences]





# final test
r_values <- seq(0, 50, length.out = 70)
# basic function overflow
#Lcross_neoStroma <- spatstat.explore::Lcross(pp_object, i = "neoplastic", j = "stromal", correction = "border", r = r_values)

# Downsampling random
# Lcross_neoStroma <- Lcross_ds(pp_object, mark_i = "neoplastic", mark_j = "stromal", 
#                               downsample = TRUE, n_down = 40000,
#                               correction = "border", r = r_values)

# Sketch downsampling
Lcross_neoStroma <- Lcross_sketch(
  pp_object2, 
  mark_i = "neoplastic", 
  mark_j = "stromal", 
  # mark_column = "merged_type", 
  downsample = TRUE, 
  n_threads = 9,
  n_down = 10000,
  correction = "border", 
  r = r_values
)


# 2 processi

keep_i <- sketch_spatial_R(pp_object2,   which(marks(pp_object2) == "neoplastic"), 10000, 8)
keep_j <- sketch_spatial_R(pp_object2, which(marks(pp_object2) == "stromal"), 10000, 8)

# Subset the ppp object
ppp_sub <- pp_object2[c(keep_i, keep_j)]

plot(ppp_sub)

ppp_sub

# second
keep_1 <- sketch_spatial_R(pp_object, which(marks(pp_object) == "neoplastic"), 8000, 8)
keep_2 <- sketch_spatial_R(pp_object, which(marks(pp_object) == "stromal"), 8000, 8)

# Subset the ppp object
ppp_sub2 <- pp_object[c(keep_1, keep_2)]
plot(ppp_sub2)
head(pp_object$y)

# -----

# con mbkmeans

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("mbkmeans")
library(mbkmeans)

ind_i <- which(marks(pp_object2) == "neoplastic")
ind_j <- which(marks(pp_object2) == "stromal")
coords <- cbind(x = pp_object2$x[ind_i], y = pp_object2$y[ind_i])

bsize <- blocksize(t(coords))
bsize

mbkm <- mbkmeans(t(coords), clusters = 10000, batch_size = bsize)


# ---- kernel est

pp_neo <- subset(pp_object2, marks == "neoplastic")
pp_neo$marks <- factor(pp_neo$marks)
levels(pp_neo$marks)

# at = "points" stima leaveoneout con somma pari a uno per ogni punto del ppp
# at = "pixels" in base alla risoluzione implementata dalle funzioni ogni punto viene assegnato al pixel più vicino
# di fatto creando una griglia in base al valore dell'intensità stimata con un kernel "grossolano" - quello che mi serve
kernint <- density.ppp(pp_neo, at = "pixels")
dim(kernint)
plot(kernint)

# campioniamo con questa cazzo! 

pints <- kernint[pp_neo]
hist(pints)
summary(pints)

# 3. Define how many points you want to keep
n_target <- 10000

# 4. Sample the indices of the points
# prob = point_intensities uses the local density as the sampling weight
sampled_indices <- sample(
  x = seq_len(npoints(pp_neo)), # Numbers from 1 to total number of points
  size = n_target,         # How many points to keep
  replace = FALSE,         # Don't pick the same point twice
  prob = pints # Probability weights based on KDE
)

# with the inverse
sampled_indices2 <- sample(
  x = seq_len(npoints(pp_neo)), # Numbers from 1 to total number of points
  size = n_target,         # How many points to keep
  replace = FALSE,         # Don't pick the same point twice
  prob = 1/pints # Probability weights based on KDE
)

# 5. Create the new downsampled point pattern
pp_sub <- pp_neo[sampled_indices]
pp_sub2 <- pp_neo[sampled_indices2]


# 1. Convert spatstat 'ppp' objects to data frames
df_original <- as.data.frame(pp_neo)
df_downsampled <- as.data.frame(pp_sub)
df_downsampled2 <- as.data.frame(pp_sub2)

# 2. Plot with ggplot2
plots <- lapply(list(df_downsampled, df_downsampled2), function(df) {
  ggplot() +
  # Plot original points in light grey in the background
  geom_point(data = df_original, aes(x = x, y = y), 
             color = "grey80", size = 0.5, alpha = 0.5) +
  
  # Plot the downsampled points in red on top
  geom_point(data = df, aes(x = x, y = y), 
             color = "red", size = 0.5, alpha = 0.8) +
  
  # Use coord_fixed() so the spatial aspect ratio isn't distorted!
  coord_fixed() +
  
  # Clean up the theme
  theme_minimal() +
  labs(title = "Spatial Downsampling (Inverse KDE)",
       subtitle = "Grey: Original | Red: Downsampled Sketch",
       x = "X", 
       y = "Y")
})
plots[[1]] + plots[[2]]


# test check ---- 

file_paths <- list.files(path = "/projects/shared/TCGA_data/results/Lfun_KDE", 
                         pattern = "\\.rds$",  # Only match .csv files
                         full.names = TRUE)    # Get the full file path!
lest_test <- imap(file_paths, \(file, name){
  name <- readRDS(file)
})

# 1. Convert the list of 'fv' objects to a single dataframe
# .id = "sample" creates a column to identify each line for coloring
df_lest <- bind_rows(lapply(lest_test, as.data.frame), .id = "sample")

# 2. Plot using ggplot2
ggplot(df_lest, aes(x = r, y = border - r, color = sample)) +
  # Dropped alpha to 0.1 so you can actually see the density of 2,000 lines
  geom_line(alpha = 0.8, size = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  guides(color = "none") +   # <--- The absolute bulletproof way to remove the color legend
  labs(title = "Centered L-function (2000 Samples)",
       x = "Distance (r)",
       y = expression(L(r) - r))


# pcf estimation

pcf_test <- map(lest_test, \(lest){
  pcf.fv(lest, method = "b")
})

df_pcf <- bind_rows(lapply(pcf_test, as.data.frame), .id = "sample")
head(df_pcf)
# 2. Plot using ggplot2
ggplot(df_pcf, aes(x = r, y = pcf, color = sample)) +
  # Dropped alpha to 0.1 so you can actually see the density of 2,000 lines
  geom_line(alpha = 0.8, size = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  guides(color = "none") +   # <--- The absolute bulletproof way to remove the color legend
  labs(title = "Centered L-function (2000 Samples)",
       x = "Distance (r)",
       y = expression(L(r) - r))
