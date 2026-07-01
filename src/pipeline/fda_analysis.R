suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(fda))
source("~/ppp_TCGA/src/plots/plot_funs.R")
source("~/ppp_TCGA/src/utils/ppp_fun.R")

# ==========================================
# 2. Define Paths
# ==========================================
input_file <- "/projects/shared/TCGA_data/results/combined_Lcross_KDE.csv"


# ==========================================
# 3. Read Data and Apply Function
# ==========================================
cat("Reading data from:", input_file, "\n")
dt <- fread(input_file)
cat("Successfully loaded data with", nrow(dt), "rows and", ncol(dt), "columns.\n")

cat("Applying fsmooth function...\n")
smoothed_Lfun <- fsmooth(dt, M = 6, genLfun.fd = 4, centrata = FALSE)

# fPCA 
#' number of basis function: length(r) = 70,  M = 6
#' K = 70 + 6 − 2 = 74 
fpca <- pca.fd(smoothed_Lfun$fd, nharm = 70, smoothed_Lfun$fdPar)
rownames(fpca$scores) <- colnames(dt)[-1] # sample names without r_values

# clustering
set.seed(123)
fpca_km <- kmeans(fpca$scores[,c(1,60)], 6)
fpca_hc <- hclust(dist(fpca$scores[,c(1,60)]), method = "ward.D2")
memb_hc <- cutree(fpca_hc, k = 6)

table(fpca_km$cluster)
table(memb_hc)

# plot of the L funs
plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun, 
  clusters_vec = fpca_km$cluster, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored


# PC scores
par(mfrow=c(1,2))
plot(fpca$scores[,1], fpca$scores[,2], col = fpca_km$cluster, pch = 3, main = "kmeans")
abline(v = 0, lty = "dashed", col = "grey")
abline(h = 0, lty = "dashed", col = "grey")
plot(fpca$scores[,1], fpca$scores[,2], col = memb_hc, pch = 3, main = "Ward.d2")
abline(v = 0, lty = "dashed", col = "grey")
abline(h = 0, lty = "dashed", col = "grey")
par(mfrow=c(1,1))

cat("Done!\n")
