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

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait        <- readRDS("output/TPDs_fish.rds")
pca_trait        <- readRDS("output/pca_trait.rds")
MatriceFish      <- read.csv("output/MatriceFish.csv")
TraitFishImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# Remove the artificial “all” row
MatriceFish <- MatriceFish[rownames(MatriceFish) != "all", , drop = FALSE]

# ------------------------------------------------------------------------------
# Threat categories
# ------------------------------------------------------------------------------

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
# Null model: FRic reduction by threat category
# ------------------------------------------------------------------------------

res_FRic_threat <- calc_FRic_by_threat( # long, recommended skip
  MatriceFish,
  TPDsp    = tpd_trait,
  threatsp = threatsp,
  nrep     = 999
)

# saveRDS(res_FRic_threat, "output/res_FRic_threat.rds")
res_FRic_threat <- readRDS("output/res_FRic_threat.rds")

# ------------------------------------------------------------------------------
# SES computation
# ------------------------------------------------------------------------------

res_SES <- calc_SES_table(res_FRic_threat) # Table S1 b

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(res_SES, "output/res_FRic_threat_SES.rds")
res_FRic_threat_SES <- readRDS("output/res_FRic_threat_SES.rds")
