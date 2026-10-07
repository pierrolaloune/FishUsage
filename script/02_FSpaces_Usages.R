# --- 02_FSpaces_Usages --------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")

pca_trait$pca_object$scores[, 1]   <- -pca_trait$pca_object$scores[, 1]
pca_trait$pca_object$loadings[, 1] <- -pca_trait$pca_object$loadings[, 1]

usage_cols      <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "All uses")
pc_combinations <- list(PC1PC2 = c(1, 2), PC3PC4 = c(3, 4))

# --- Global functional spaces -------------------------------------------------

FS_global_PC1PC2 <- funspace(x = pca_trait$pca_object, PCs = c(1, 2), n_divisions = 300, threshold = 0.999)
summary(FS_global_PC1PC2)

FS_global_PC3PC4 <- funspace(x = pca_trait$pca_object, PCs = c(3, 4), n_divisions = 300, threshold = 0.999)
summary(FS_global_PC3PC4)

# --- Functional spaces [LONG] -------------------------------------------------

# funspace_results <- list()
#
# for (use in usage_cols) {
#   for (pc_name in names(pc_combinations)) {
#     funspace_results[[paste0("FS_", gsub(" ", "", use), "_", pc_name)]] <- funspace(
#       x           = pca_trait$pca_object,
#       group.vec   = as.factor(pca_trait$uses[[use]]),
#       PCs         = pc_combinations[[pc_name]],
#       n_divisions = 300,
#       threshold   = 0.999
#     )
#   }
# }
#
# saveRDS(funspace_results, "output/funspace_results.rds")

funspace_results <- readRDS("output/funspace_results.rds")

# --- Session ------------------------------------------------------------------

sessionInfo()
