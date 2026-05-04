# ------------------------------------------------------------------------------
# Script : 03_PCA_mean_trait_value
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads fish TPDs, PCA outputs, and the community-by-use matrix, then
# runs a null model to simulate mean PCA scores per human-use community. It
# computes SES for observed mean PCA values against simulated distributions and
# prepares long-format data for visualization. It produces (i) histograms of
# simulated PCA distributions with observed values overlaid, (ii) a heatmap of
# signed trait loadings across PCA axes, and (iii) a PCA correlation circle for
# PC1–PC2 based on loadings and explained variance.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait     <- readRDS("output/TPDs_fish.rds")
pca_trait     <- readRDS("output/PCA_fish.rds")
pca_trait_use <- readRDS("output/pca_trait.rds")

MatriceFish <- read.csv("output/MatriceFish.csv")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# ------------------------------------------------------------------------------
# Null model of PCA mean trait value
# ------------------------------------------------------------------------------

resultats_null <- generate_null_means( # long recommended skip
  pca_trait      = pca_trait,
  MatriceFish    = MatriceFish,
  nb_simulations = 999
)

resultats_null$all <- NULL

# saveRDS(resultats_null, "output/PCA_mean_trait_values_results.rds")
PCA_mean_trait_values_results <- readRDS("output/PCA_mean_trait_values_results.rds")
# Invert PC1 so that higher values indicate larger body size (biological interpretability)
PCA_mean_trait_values_results <- lapply(PCA_mean_trait_values_results, function(use) {
  use$observed["Comp.1"]    <- -use$observed["Comp.1"]
  use$simulated[, "Comp.1"] <- -use$simulated[, "Comp.1"]
  use
})

# ------------------------------------------------------------------------------
# SES computation
# ------------------------------------------------------------------------------

PCA_mean_trait_values_SES <- get_SES_from_PCA_results( # Table S1
  results_list = PCA_mean_trait_values_results
)

# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")
