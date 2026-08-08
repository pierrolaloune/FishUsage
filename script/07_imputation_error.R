# ------------------------------------------------------------------------------
# Script : 07_imputation_error
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads raw (missing) and imputed trait tables, the PCA/use object,
# scaling parameters, and a phylogeny. It then evaluates phylogeny-informed
# missForest imputation error using bootstrap masking simulations, summarizing
# NRMSE across iterations and saving/loading the resulting output object.
#
# Principle: species with complete traits are masked using a real missingness
# pattern, re-imputed, and projected back into the reference PCA. The distance
# between the projected and the true position gives the error, expressed as a
# percentage of the range of each axis.
#
# The evaluation itself is commented out; its result is stored in output/ and
# reloaded at the end, so the script runs end to end as it is.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
traitsData <- read.table("dataPrepared/Fish/TraitFishMissing.txt") %>%
  dplyr::select(-IUCN)

traitsDataImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt") %>%
  dplyr::select(-IUCN)

selectedTraits <- colnames(traitsDataImputed)

pca_trait <- readRDS("output/pca_trait.rds")
PCAmodel  <- pca_trait$pca_object
traitPCA  <- pca_trait$traits_scores[, 1:4]

# ---- Centring and scaling used by the reference PCA ----
meanInputed <- attr(pca_trait$traits_scaled, "scaled:center")
sdInputed   <- attr(pca_trait$traits_scaled, "scaled:scale")

phylogeny <- readRDS("dataOriginal/FishMORPH_Phylogeny.rds")

# ------------------------------------------------------------------------------
# Imputation error settings
# ------------------------------------------------------------------------------

nboot_val   <- 100    # number of bootstrap iterations
perc_val    <- 0.1    # proportion of species masked each iteration
npcoa_val   <- 2      # number of phylogenetic PCoA axes
seed_val    <- 123
ref_max_val <- 1000   # maximum complete species used as reference
ntree_val   <- 30     # trees per random forest
maxiter_val <- 2      # missForest iterations

# ------------------------------------------------------------------------------
# Run imputation error model  [LONG]
# ------------------------------------------------------------------------------

# LONG: one full missForest run per bootstrap iteration (100 by default).

# cat("\n=== STARTING PHYLOGENETIC IMPUTATION ERROR ===\n")
# set.seed(seed_val)
#
# res <- evaluate_imputation_phylo(
#   traitsData        = traitsData,
#   traitsDataImputed = traitsDataImputed,
#   selectedTraits    = selectedTraits,
#   meanInputed       = meanInputed,
#   sdInputed         = sdInputed,
#   traitPCA          = traitPCA,
#   PCAmodel          = PCAmodel,
#   phylogeny         = phylogeny,
#   dimensions        = 1:4,
#   percImpute        = perc_val,
#   nboot             = nboot_val,
#   npcoa             = npcoa_val,
#   ncores            = 1,              # force sequential mode
#   ref_complete_max  = ref_max_val,
#   ntree             = ntree_val,
#   maxiter           = maxiter_val,
#   seed              = seed_val
# )
#
# cat("\n=== ANALYSIS COMPLETE ===\n")
# print(res$summary)
#
# saveRDS(res, "output/NRMSE_results.rds")

# ------------------------------------------------------------------------------
# Results
# ------------------------------------------------------------------------------

# ---- Recommended: load the saved result ----
NMRSE_results <- readRDS("output/NRMSE_results.rds")
NMRSE_summary <- NMRSE_results$summary
