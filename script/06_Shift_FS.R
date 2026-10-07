# --- 06_Shift_FS --------------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")

pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

# --- TPD [LONG] ---------------------------------------------------------------

sd_traits <- sqrt(diag(ks::Hpi.diag(pca_trait$traits_scores[, c(1, 2)])))

# TPD_2D <- TPD::TPDsMean(
#   species      = rownames(pca_trait$traits_scores),
#   means        = pca_trait$traits_scores[, c(1, 2)],
#   sds          = matrix(rep(sd_traits, nrow(pca_trait$traits_scores)), byrow = TRUE, ncol = 2),
#   covar        = FALSE,
#   alpha        = 0.95,
#   samples      = NULL,
#   trait_ranges = NULL,
#   n_divisions  = 200,
#   tolerance    = 0.05
# )
#
# saveRDS(TPD_2D, "output/TPD_2D.rds")

TPD_2D <- readRDS("output/TPD_2D.rds")

# --- Session ------------------------------------------------------------------

sessionInfo()
