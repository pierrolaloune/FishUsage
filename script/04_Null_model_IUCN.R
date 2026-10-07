# --- 04_Null_model_IUCN -------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

tpd_trait        <- readRDS("output/TPDs_fish.rds")
TraitFishImputed <- read.table("dataPrepared/Fish/TraitFishImputed.txt")

MatriceFish <- read.csv("output/MatriceFish.csv")
colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

MatriceFish <- MatriceFish[rownames(MatriceFish) != "all", , drop = FALSE]

IUCN_levels <- list(
  CR             = c("CR"),
  CR_EN          = c("CR", "EN"),
  CR_EN_VU       = c("CR", "EN", "VU"),
  CR_EN_VU_NT    = c("CR", "EN", "VU", "NT"),
  CR_EN_VU_NT_DD = c("CR", "EN", "VU", "NT", "DD")
)

threatsp <- lapply(IUCN_levels, function(levels) {
  rownames(TraitFishImputed)[TraitFishImputed$IUCN %in% levels]
})

# --- Null model [LONG] --------------------------------------------------------

# res_FRic_threat <- calc_FRic_by_threat(MatriceFish, TPDsp = tpd_trait, threatsp = threatsp, nrep = 999)
# saveRDS(res_FRic_threat, "output/res_FRic_threat.rds")

res_FRic_threat <- readRDS("output/res_FRic_threat.rds")

# --- SES ----------------------------------------------------------------------

res_SES <- calc_SES_table(res_FRic_threat)
print(res_SES)

# --- Save ---------------------------------------------------------------------

# saveRDS(res_SES, "output/res_FRic_threat_SES.rds")

# --- Session ------------------------------------------------------------------

sessionInfo()
