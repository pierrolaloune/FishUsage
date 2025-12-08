# ============================================================
# Script: 06_Shift_FS.R
# Purpose: Quantify functional space shifts after removing
#          threatened species, using 2D TPD and usage categories
# ============================================================

################################################################################
# DATA PREPARATION
################################################################################

pca_trait <- readRDS("output/pca_trait.rds")

# Flip PC1 for readability
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

# SD estimates for 2D kernel smoothing
sd_traits <- sqrt(diag(
  ks::Hpi.diag(pca_trait$traits_scores[, c(1, 2)])
))

################################################################################
# 2D TPD CALCULATION (LONG)
################################################################################

TPD_2D <- TPDsMean( # long recommended skip
  species = rownames(pca_trait$traits_scores),
  means   = pca_trait$traits_scores[, c(1, 2)],
  sds     = matrix(rep(sd_traits, nrow(pca_trait$traits_scores)),
                   byrow = TRUE, ncol = 2),
  covar        = FALSE,
  alpha        = 0.95,
  samples      = NULL,
  trait_ranges = NULL,
  n_divisions  = 200,
  tolerance    = 0.05
)

# saveRDS(TPD_2D, "output/TPD_2D.rds")

################################################################################
# LOAD OBJECTS
################################################################################

TPDs_fish <- readRDS("output/TPD_2D.rds")
uni_clean  <- readRDS("output/uni_clean.rds")

IUCN <- uni_clean %>%
  dplyr::select(species, IUCN)

################################################################################
# PLOT PARAMETERS
################################################################################

limX   <- c(-7, 7)
limY   <- c(-7, 7)
usages <- c("Fisheries", "Aquaculture", "Aquarium",
            "Game fish", "Bait", "All uses")

################################################################################
# FS SHIFT PLOTS (PER USE)
################################################################################

for (u in usages) {
  plot_functional_shift_by_usage(usage_name = u)
}

################################################################################
# LEGEND (COLORBAR)
################################################################################

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
  zlim        = c(Min, Max),
  legend.only = TRUE,
  col         = rampcols,
  breaks      = rampbreaks,
  horizontal  = FALSE,
  legend.width = 1.2,
  legend.mar   = 4,
  axis.args    = list(
    at     = seq(-0.3, 0.3, by = 0.1),
    labels = paste0(seq(-30, 30, by = 10), "%")
  )
)

dev.off()
