# ============================================================
# Script: 04_Null_model_IUCN.R
# Purpose: Evaluate the impact of removing threatened species
#          on functional richness (FRic) through a null model,
#          and compute corresponding SES values
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

tpd_trait         <- readRDS("output/TPDs_fish.rds")
pca_trait         <- readRDS("output/pca_trait.rds")
MatriceFish       <- read.csv("output/MatriceFish.csv")
TraitFishImputed  <- read.table("dataPrepared/Fish/TraitFishImputed.txt")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# Remove the artificial “all” row
MatriceFish <- MatriceFish[rownames(MatriceFish) != "all", , drop = FALSE]

################################################################################
# THREAT CATEGORIES
################################################################################

species_names <- rownames(TraitFishImputed)

IUCN_levels <- list(
  CR              = c("CR"),
  CR_EN           = c("CR", "EN"),
  CR_EN_VU        = c("CR", "EN", "VU"),
  CR_EN_VU_NT     = c("CR", "EN", "VU", "NT"),
  CR_EN_VU_NT_DD  = c("CR", "EN", "VU", "NT", "DD")
)

threatsp <- lapply(IUCN_levels, function(levels) {
  species_names[TraitFishImputed$IUCN %in% levels]
})

################################################################################
# NULL MODEL: FRIC REDUCTION BY THREAT CATEGORY
################################################################################

res_FRic_threat <- calc_FRic_by_threat( # long, recommended skip
  MatriceFish,
  TPDsp     = tpd_trait,
  threatsp  = threatsp,
  nrep      = 999
)

# saveRDS(res_FRic_threat, "output/res_FRic_threat.rds")
res_FRic_threat <- readRDS("output/res_FRic_threat.rds")

################################################################################
# SES COMPUTATION
################################################################################

res_SES <- calc_SES_table(res_FRic_threat) # Table S1 b

# saveRDS(res_SES, "output/res_FRic_threat_SES.rds")
res_FRic_threat_SES <- readRDS("output/res_FRic_threat_SES.rds")
