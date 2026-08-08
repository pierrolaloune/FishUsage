# ------------------------------------------------------------------------------
# Script : 04_Null_model_IUCN
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads fish TPDs, the PCA/use object, the community-by-use matrix,
# and the imputed trait table with IUCN categories. It defines nested threat
# categories (from CR alone up to CR+EN+VU+NT+DD), builds corresponding species
# lists, and runs a null model to quantify the reduction in functional richness
# (FRic) when threatened species are removed. It then computes SES values from
# the null distributions and loads saved outputs for downstream reporting.
#
# Output: Table S1b (SES of the FRic loss per use x threat category).
#
# The null model is commented out; its result is stored in output/ and reloaded
# immediately afterwards, so the script runs end to end as it is.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait        <- readRDS("output/TPDs_fish.rds")
pca_trait        <- readRDS("output/pca_trait.rds")
MatriceFish      <- read.csv("output/MatriceFish.csv")
TraitFishImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt")

# ---- Restore the species names and the row names lost by read.csv ----
colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# ---- Drop the artificial "all" row: only real uses are tested here ----
MatriceFish <- MatriceFish[rownames(MatriceFish) != "all", , drop = FALSE]

# ------------------------------------------------------------------------------
# Threat categories
# ------------------------------------------------------------------------------

# Nested scenarios: each level adds one IUCN category to the previous one, which
# simulates an increasingly severe loss of threatened species.

species_names <- rownames(TraitFishImputed)

IUCN_levels <- list(
  CR             = c("CR"),
  CR_EN          = c("CR", "EN"),
  CR_EN_VU       = c("CR", "EN", "VU"),
  CR_EN_VU_NT    = c("CR", "EN", "VU", "NT"),
  CR_EN_VU_NT_DD = c("CR", "EN", "VU", "NT", "DD")
)

threatsp <- lapply(IUCN_levels, function(levels) {
  species_names[TraitFishImputed$IUCN %in% levels]
})

# ------------------------------------------------------------------------------
# Null model: FRic reduction by threat category  [LONG]
# ------------------------------------------------------------------------------

# LONG: 999 random draws for each of the 5 uses x 5 threat categories, each draw
# rebuilding a community TPD.

# res_FRic_threat <- calc_FRic_by_threat(
#   MatriceFish,
#   TPDsp    = tpd_trait,
#   threatsp = threatsp,
#   nrep     = 999
# )
#
# saveRDS(res_FRic_threat, "output/res_FRic_threat.rds")

# ---- Recommended: load the saved result ----
res_FRic_threat <- readRDS("output/res_FRic_threat.rds")

# ------------------------------------------------------------------------------
# SES computation
# ------------------------------------------------------------------------------

# ---- Table S1b ----
res_SES <- calc_SES_table(res_FRic_threat)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(res_SES, "output/res_FRic_threat_SES.rds")
res_FRic_threat_SES <- readRDS("output/res_FRic_threat_SES.rds")
