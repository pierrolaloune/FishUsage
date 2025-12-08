# ============================================================
# Script: 02_FSpaces_Usages.R
# Purpose: Build functional spaces (FS) for the full fish dataset
#          and for each human use category using PCA + funspace
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

pca_trait <- readRDS("output/pca_trait.rds")
tpd_trait <- readRDS("output/TPDs_fish.rds")

################################################################################
# GLOBAL FUNCTIONAL SPACES (FULL COMMUNITY)
################################################################################

# Flip PC1 for interpretability
pca_trait$pca_object$scores[, 1]    <- -pca_trait$pca_object$scores[, 1]
pca_trait$pca_object$loadings[, 1]  <- -pca_trait$pca_object$loadings[, 1]

# Global FS: PC1 + PC2
FS_global_PC1PC2 <- funspace(
  x           = pca_trait$pca_object,
  PCs         = c(1, 2),
  n_divisions = 300,
  threshold   = 0.999
)
summary(FS_global_PC1PC2)

# Global FS: PC3 + PC4
FS_global_PC3PC4 <- funspace(
  x           = pca_trait$pca_object,
  PCs         = c(3, 4),
  n_divisions = 300,
  threshold   = 0.999
)
summary(FS_global_PC3PC4)

################################################################################
# FUNCTIONAL SPACES BY HUMAN USE CATEGORY
################################################################################

usage_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "Bait", "All uses")

pc_combinations <- list(
  PC1PC2 = c(1, 2),
  PC3PC4 = c(3, 4)
)

funspace_results <- list()

# Long computation
for (use in usage_cols) {
  use_clean <- gsub(" ", "", use)
  
  for (pc_name in names(pc_combinations)) {
    pcs <- pc_combinations[[pc_name]]
    
    res <- funspace(
      x           = pca_trait$pca_object,
      group.vec   = as.factor(pca_trait$uses[[use]]),
      PCs         = pcs,
      n_divisions = 300,
      threshold   = 0.999
    )
    
    result_name <- paste0("FS_", use_clean, "_", pc_name)
    funspace_results[[result_name]] <- res
  }
}

# saveRDS(funspace_results, "output/funspace_results.rds")
# saveRDS(funspace_results, "output/funspace_results_compressed.rds", compress = "xz")
