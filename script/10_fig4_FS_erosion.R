# ------------------------------------------------------------------------------
# Script : 10_fig4_FS_erosion.R
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads the PCA/use object, flips PC1 for readability, estimates
# kernel bandwidths for 2D trait-space smoothing, and computes a 2D TPD across
# species (optionally, long computation). It then loads precomputed 2D TPD and
# cleaned species/IUCN data, iterates over human-use categories to plot functional
# space shifts after removing threatened species, and finally generates a shared
# diverging colorbar legend for the shift maps.
#
# Outputs: figures/FS_shift_<usage>.jpg and figures/FS_shift_legend.jpg.
# The maps themselves are drawn by plot_functional_shift_by_usage(), defined in
# 000_functions.R. That function reads pca_trait, IUCN, TPDs_fish, limX and limY
# from the global environment, so all of them must exist before the loop below.

# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")

# ---- Flip PC1 so that large-bodied species sit on the right ----
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

# ---- Bandwidths for the 2D kernel smoothing (plug-in estimator) ----
sd_traits <- sqrt(diag(
  ks::Hpi.diag(pca_trait$traits_scores[, c(1, 2)])
))

# ------------------------------------------------------------------------------
# 2D TPD calculation  [LONG]
# ------------------------------------------------------------------------------

# LONG: one bivariate density per species on a 200 x 200 grid.

# TPD_2D <- TPDsMean(
#   species = rownames(pca_trait$traits_scores),
#   means   = pca_trait$traits_scores[, c(1, 2)],
#   sds     = matrix(
#     rep(sd_traits, nrow(pca_trait$traits_scores)),
#     byrow = TRUE, ncol = 2
#   ),
#   covar        = FALSE,
#   alpha        = 0.95,
#   samples      = NULL,
#   trait_ranges = NULL,
#   n_divisions  = 200,
#   tolerance    = 0.05
# )
#
# saveRDS(TPD_2D, "output/TPD_2D.rds")

# ------------------------------------------------------------------------------
# Load objects
# ------------------------------------------------------------------------------

# ---- Recommended: load the saved 2D TPD ----
TPDs_fish <- readRDS("output/TPD_2D.rds")
uni_clean <- readRDS("output/uni_clean.rds")

IUCN <- uni_clean %>%
  dplyr::select(species, IUCN)

# ------------------------------------------------------------------------------
# Plot parameters
# ------------------------------------------------------------------------------

# Axis limits are fixed so that every map shares the same frame.
limX   <- c(-7, 7)
limY   <- c(-7, 7)
usages <- c("Fisheries", "Aquaculture", "Aquarium",
            "Game fish", "Bait", "All uses")

# ------------------------------------------------------------------------------
# FS shift plots (per use)
# ------------------------------------------------------------------------------

for (u in usages) {
  plot_functional_shift_by_usage(usage_name = u)
}

# ------------------------------------------------------------------------------
# Legend (colorbar)
# ------------------------------------------------------------------------------

# Diverging scale centred on 0, drawn once and shared by every map above.

Min    <- -0.3
Max    <- 0.3
Thresh <- 0
ncol   <- 1000

ColorRamp <- rev(scico::scico(n = ncol, palette = "vik"))

nHalf <- 500

rc1 <- grDevices::colorRampPalette(ColorRamp[1:nHalf], space = "Lab")(nHalf)
rc2 <- grDevices::colorRampPalette(ColorRamp[(nHalf + 1):ncol], space = "Lab")(nHalf)
rampcols <- c(rc1, rc2)

rampbreaks <- c(
  seq(Min, Thresh, length.out = nHalf + 1),
  seq(Thresh, Max, length.out = nHalf + 1)[-1]
)

jpeg("figures/FS_shift_legend.jpg", width = 500, height = 1600, res = 300)

par(mar = c(4, 5, 2, 2))

fields::image.plot(
  zlim         = c(Min, Max),
  legend.only  = TRUE,
  col          = rampcols,
  breaks       = rampbreaks,
  horizontal   = FALSE,
  legend.width = 1.2,
  legend.mar   = 4,
  axis.args    = list(
    at     = seq(-0.3, 0.3, by = 0.1),
    labels = paste0(seq(-30, 30, by = 10), "%")
  )
)

dev.off()

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(TPD_2D, "output/TPD_2D.rds")
