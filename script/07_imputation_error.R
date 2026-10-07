# --- 07_imputation_error ------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

traitsData        <- read.table("dataPrepared/Fish/TraitFishMissing.txt") %>% dplyr::select(-IUCN)
traitsDataImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt") %>% dplyr::select(-IUCN)
selectedTraits    <- colnames(traitsDataImputed)

pca_trait <- readRDS("output/pca_trait.rds")
PCAmodel  <- pca_trait$pca_object
traitPCA  <- pca_trait$traits_scores[, 1:4]

meanInputed <- attr(pca_trait$traits_scaled, "scaled:center")
sdInputed   <- attr(pca_trait$traits_scaled, "scaled:scale")

phylogeny <- readRDS("dataPrepared/Fish/FishMORPH_Phylogeny.rds")

nboot_val   <- 100
perc_val    <- 0.1
npcoa_val   <- 2
seed_val    <- 123
ref_max_val <- 1000
ntree_val   <- 30
maxiter_val <- 2

# --- Imputation error [LONG] --------------------------------------------------

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
#   ncores            = 1,
#   ref_complete_max  = ref_max_val,
#   ntree             = ntree_val,
#   maxiter           = maxiter_val,
#   seed              = seed_val
# )
#
# saveRDS(res, "output/NRMSE_results.rds")

NMRSE_results <- readRDS("output/NRMSE_results.rds")
NMRSE_summary <- NMRSE_results$summary
print(NMRSE_summary)

# --- Session ------------------------------------------------------------------

sessionInfo()
