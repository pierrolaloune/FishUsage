# ============================================================
# Script: 07_imputation_error.R
# Purpose: Evaluate imputation error using phylogeny-informed
#          missForest and bootstrap simulations (NRMSE)
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

traitsData <- read.table("dataPrepared/Fish/TraitFishMissing.txt") %>%
  dplyr::select(-IUCN)

traitsDataImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt") %>%
  dplyr::select(-IUCN)

selectedTraits <- colnames(traitsDataImputed)

pca_trait <- readRDS("output/pca_trait.rds")
PCAmodel  <- pca_trait$pca_object
traitPCA  <- pca_trait$traits_scores[, 1:4]

meanInputed <- attr(pca_trait$traits_scaled, "scaled:center")
sdInputed   <- attr(pca_trait$traits_scaled, "scaled:scale")

phylogeny <- readRDS("dataOriginal/FishMORPH_Phylogeny.rds")

################################################################################
# IMPUTATION ERROR SETTINGS
################################################################################

nboot_val   <- 30     # number of bootstrap iterations
perc_val    <- 0.1    # proportion of species masked each iteration
npcoa_val   <- 2      # number of phylogenetic PCoA axes
seed_val    <- 123
ref_max_val <- 1000   # maximum complete species used as reference
ntree_val   <- 30
maxiter_val <- 2

################################################################################
# RUN IMPUTATION ERROR MODEL
################################################################################

cat("\n=== STARTING PHYLOGENETIC IMPUTATION ERROR ===\n")
set.seed(seed_val)

res <- evaluate_imputation_phylo(
  traitsData        = traitsData,
  traitsDataImputed = traitsDataImputed,
  selectedTraits    = selectedTraits,
  meanInputed       = meanInputed,
  sdInputed         = sdInputed,
  traitPCA          = traitPCA,
  PCAmodel          = PCAmodel,
  phylogeny         = phylogeny,
  dimensions        = 1:4,
  percImpute        = perc_val,
  nboot             = nboot_val,
  npcoa             = npcoa_val,
  ncores            = 1,              # force sequential mode
  ref_complete_max  = ref_max_val,
  ntree             = ntree_val,
  maxiter           = maxiter_val,
  seed              = seed_val
)

cat("\n=== ANALYSIS COMPLETE ===\n")
print(res$summary)

# saveRDS(res, "output/NRMSE_results.rds")
NMRSE_results <- readRDS("output/NRMSE_results.rds")
