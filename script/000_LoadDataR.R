# ============================================================
# Script: 000_LoadDataR.R
# Purpose: Load and prepare trait, phylogeny, IUCN and use data
#          for freshwater fishes, and build matrices for TPD
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

# Raw trait and phylogeny from FishMORPH
trait     <- readRDS("dataOriginal/FishMORPH_Traits.rds")
phylogeny <- readRDS("dataOriginal/FishMORPH_Phylogeny.rds")
traitNames <- colnames(trait)[-c(1:6)]

################################################################################
# TRAITS + FISHBASE LENGTH/WEIGHT
################################################################################

# Species list
list_sp <- gsub("\\.", " ", as.character(trait$Genus.species))

# Length-weight from FishBase
lgtwgt <- as.data.frame(length_weight(list_sp))
lgtwgt <- lgtwgt[!is.na(lgtwgt$a) & !is.na(lgtwgt$b), ]

# Best a, b coefficients per species
ab <- lgtwgt %>%
  dplyr::group_by(Species) %>%
  dplyr::slice_max(order_by = CoeffDetermination, with_ties = FALSE, na_rm = TRUE) %>%
  dplyr::select(Species, a, b) %>%
  dplyr::distinct()

# Additional traits from FishBase
speciesInfo <- species(list_sp) %>% data.table::data.table()
speciesInfoSub <- speciesInfo[, .(Species, Fresh, LongevityWild, Length, Weight, LTypeMaxM)]
speciesInfoSub <- unique(speciesInfoSub)
speciesInfoSub$Species <- gsub(" ", ".", speciesInfoSub$Species)

speciesInfoSub <- speciesInfoSub %>%
  dplyr::mutate(
    Length = ifelse(LTypeMaxM != "SL", NA, Length),
    Weight2 = ifelse(
      Species %in% gsub(" ", ".", rownames(ab)),
      ab[gsub("\\.", " ", Species), "a"] * Length^ab[gsub("\\.", " ", Species), "b"],
      NA
    )
  ) %>%
  dplyr::mutate(
    Weight = ifelse(is.na(Weight) & !is.na(Weight2), Weight2, Weight)
  )

# Merge FishMORPH traits with FishBase length/weight
fishTraits <- merge(
  trait[, 6:15],
  speciesInfoSub[, .(Species, Length, Weight)],
  by.x = "Genus.species",
  by.y = "Species",
  all.x = TRUE
)

fishTraits <- fishTraits %>%
  dplyr::mutate(dplyr::across(-Genus.species, ~ log10(.x + 1))) %>%
  dplyr::rename(species = Genus.species) %>%
  dplyr::mutate(species = gsub("\\.", "_", species))

################################################################################
# FILTER FRESHWATER SPECIES
################################################################################

spToKeep <- fishTraits %>%
  dplyr::select(species) %>%
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::pull() %>%
  rfishbase::species() %>%
  data.table::as.data.table() %>%
  dplyr::filter(Fresh == 1) %>%
  dplyr::select(Species) %>%
  dplyr::mutate(Species = gsub(" ", "_", Species)) %>%
  dplyr::pull()

fishTraits <- fishTraits %>%
  dplyr::filter(species %in% spToKeep)

# If needed:
# dir.create("dataPrepared/Fish", showWarnings = FALSE, recursive = TRUE)
# write.table(fishTraits, "dataPrepared/Fish/fishTraitsMissing.txt")

fishTraits <- read.table("dataPrepared/Fish/fishTraitsMissing.txt")

################################################################################
# PHYLOGENETIC PCoA
################################################################################

phylogeny$tip.label <- gsub("\\.", "_", phylogeny$tip.label)
phylogeny <- drop.tip(phylogeny, setdiff(phylogeny$tip.label, fishTraits$species))
phylogenyTraits <- phytools::force.ultrametric(phylogeny)
phylDiss <- sqrt(cophenetic(phylogenyTraits))

# Long computation alternative:
# pcoaPhyl <- cmdscale(phylDiss, k = 10)

# Pre-computed PCoA (recommended)
pcoaPhyl <- read.table(
  "dataPrepared/Fish/pcoaPhylogenyFish.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

rownames(pcoaPhyl) <- fishTraits$species
colnames(pcoaPhyl) <- paste0("Eigen.", 1:10)

# If you want to recompute and save:
# write.table(pcoaPhyl, "dataPrepared/Fish/pcoaPhylogenyFish.txt")

################################################################################
# TAXONOMIC STANDARDIZATION 
################################################################################

list_sp_raw <- fishTraits$species
list_sp_raw <- gsub("_", " ", list_sp_raw)
list_sp_raw <- stringr::str_squish(list_sp_raw)

gna_one_safe <- function(x) {
  tryCatch(
    taxize::gna_verifier(
      names         = x,
      data_sources  = 11,
      all_matches   = FALSE,
      capitalize    = TRUE,
      species_group = TRUE,
      output_type   = "table"
    ),
    error = function(e) {
      tibble::tibble(
        submittedName = x,
        matchedName   = NA_character_,
        matchType     = "Error",
        dataSourceId  = NA_real_
      )
    }
  )
}

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

species_plan <- list_sp_raw

verified_names <- progressr::with_progress({  # long
  p <- progressr::progressor(along = species_plan)
  purrr::map_dfr(species_plan, function(x) {
    p(message = x)
    gna_one_safe(x)
  })
})

new_names <- gsub(" ", "_", verified_names$matchedCanonicalSimple)

rownames(pcoaPhyl) <- new_names
fishTraits$species <- new_names
traitsAndPCOA      <- cbind(fishTraits, pcoaPhyl)

# write.table(traitsAndPCOA, "dataPrepared/Fish/traitsWithPCOA.txt")
traitsAndPCOA <- read.table(
  "dataPrepared/Fish/traitsWithPCOA.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

################################################################################
# IUCN DATA (2024) + SYNONYMS
################################################################################

species_list <- unique(stringr::str_trim(traitsAndPCOA$species))
iucn_statut <- read.csv(
  "dataOriginal/assessments.csv",
  header = TRUE,
  stringsAsFactors = FALSE
) # IUCN data 2024

iucn_clean <- iucn_statut %>%
  dplyr::mutate(scientificName = stringr::str_trim(scientificName))

species_list_spaces <- gsub("_", " ", species_list)

matched_species <- dplyr::inner_join(
  tibble::tibble(species = species_list_spaces),
  iucn_clean,
  by = c("species" = "scientificName")
)

species_to_update <- tibble::tibble(species = species_list_spaces) %>%
  dplyr::anti_join(iucn_clean, by = c("species" = "scientificName"))

synonyms_info <- rfishbase::synonyms(
  species_list = species_to_update$species,
  server       = "fishbase",
  version      = "latest",
  fields       = NULL
)

synonyms_info_filtered <- synonyms_info %>%
  dplyr::filter(!Status %in% c("misapplied name", "ambiguous synonym", "provisionally accepted name"))

synonyms_mapping <- synonyms_info_filtered %>%
  dplyr::select(Species, synonym) %>%
  dplyr::distinct()

iucn_clean <- iucn_clean %>%
  dplyr::mutate(
    scientificName = dplyr::if_else(
      scientificName %in% synonyms_mapping$Species,
      synonyms_mapping$synonym[match(scientificName, synonyms_mapping$Species)],
      scientificName
    )
  )

matched_species <- dplyr::inner_join(
  tibble::tibble(species = species_list_spaces),
  iucn_clean,
  by = c("species" = "scientificName")
)

species_to_update <- tibble::tibble(species = species_list_spaces) %>%
  dplyr::anti_join(iucn_clean, by = c("species" = "scientificName")) %>%
  as.data.frame()
rownames(species_to_update) <- species_to_update$species

# Manual corrections (pre-checked)
species_to_check <- read.csv(
  "dataOriginal/species_to_update_900_done.csv",
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE
)
colnames(species_to_check) <- c("scientificName", "redlistCategory")

iucn_clean <- iucn_clean[, c("scientificName", "redlistCategory")]

acronyms <- c(
  "Critically Endangered" = "CR",
  "Endangered"            = "EN",
  "Vulnerable"            = "VU",
  "Near Threatened"       = "NT",
  "Least Concern"         = "LC",
  "Data Deficient"        = "DD",
  "Extinct"               = "EX",
  "Extinct in the Wild"   = "EW",
  "Not Evaluated"         = "NE"
)

iucn_clean <- iucn_clean %>%
  dplyr::mutate(redlistCategory = dplyr::recode(redlistCategory, !!!acronyms)) %>%
  dplyr::bind_rows(species_to_check) %>%
  dplyr::filter(scientificName %in% species_list_spaces) %>%
  dplyr::group_by(scientificName) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

traitsAndPCOA$species <- gsub("_", " ", traitsAndPCOA$species)
traitsAndPCOA$IUCN <- iucn_clean$redlistCategory[
  match(traitsAndPCOA$species, iucn_clean$scientificName)
]

# write.table(traitsAndPCOA, "dataPrepared/Fish/traitsWithPCOAIUCN.txt")
traitsAndPCOAIUCN <- read.table(
  "dataPrepared/Fish/traitsWithPCOAIUCN.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

################################################################################
# COMPLETE TAXONOMY FROM GBIF
################################################################################

taxInfo <- pbapply::pblapply(traitsAndPCOAIUCN$species, function(sp) {
  tryCatch({
    traitdataform::get_gbif_taxonomy(
      sp,
      subspecies       = TRUE,
      higherrank       = TRUE,
      conf_threshold   = 80,
      resolve_synonyms = FALSE
    )[1, ]
  }, error = function(e) NULL)
}) %>%
  data.table::rbindlist(fill = TRUE)

################################################################################
# FINAL TABLE + IMPUTATION SETUP
################################################################################

fishTraitsPhylogenyIUCN <- data.table::data.table(traitsAndPCOAIUCN)

# write.table(fishTraitsPhylogenyIUCN, "dataPrepared/Fish/AllDataFish_clean.txt")
fishData <- read.table(
  "dataPrepared/Fish/traitsWithPCOAIUCN.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Columns for traits and imputation
columnsTraits     <- 2:(which(colnames(fishData) == "Eigen.1") - 1)
columnsImputation <- 2:(which(colnames(fishData) == "IUCN") - 1)

################################################################################
# MISSFOREST IMPUTATION
################################################################################

# Long computation:

# set.seed(123)
# imputed_forest <- missForest(xmis = fishData[, columnsImputation]) # long
# print(imputed_forest$OOBerror)
# 
# traits_names <- colnames(fishData)[columnsTraits]
# 
# fishData_imputed_forest <- fishData
# fishData_imputed_forest[traits_names] <- imputed_forest$ximp[traits_names]

# write.table(
#   fishData_imputed_forest,
#   "dataPrepared/Fish/fishData_imputed_forest.txt",
#   row.names = FALSE
# )

# recommended 
fishData_imputed_forest <- read.table(
  "dataPrepared/Fish/fishData_imputed_forest.txt",
  header = TRUE
)

################################################################################
# TRAIT SELECTION AND RENAMING
################################################################################

selectedTraits <- c(
  "EdHd", "EhBd", "JlHd", "MoBd", "BlBd", "HdBd",
  "PFiBd", "PFlBl", "CFdCPd", "Length", "Weight"
)

newTraitNames <- c(
  "es", "ep", "ms", "mp", "elo", "wid",
  "pp", "ps", "cs", "svl", "bm"
)

fishTraitsMissing <- fishData[, selectedTraits]
colnames(fishTraitsMissing) <- newTraitNames
rownames(fishTraitsMissing) <- fishData$species

fishTraitsImputed <- fishData_imputed_forest[, selectedTraits]
colnames(fishTraitsImputed) <- newTraitNames
rownames(fishTraitsImputed) <- fishData$species

# Add IUCN
fishTraitsMissing <- data.frame(fishTraitsMissing, IUCN = fishData$IUCN)
fishTraitsImputed <- data.frame(fishTraitsImputed, IUCN = fishData$IUCN)

# write.table(fishTraitsMissing, "dataPrepared/Fish/TraitFishMissing.txt")
# write.table(fishTraitsImputed, "dataPrepared/Fish/TraitFishImputed.txt")

################################################################################
# PCA + TPD (MORPHOLOGICAL SPACE)
################################################################################

# Long computation:
# results <- computePCAandTPDs(fishTraitsImputed[, !colnames(fishTraitsImputed) == "IUCN"])
# saveRDS(results, "output/All_fish.rds")

results   <- readRDS("output/All_fish.rds")
pca_trait <- results$PCA
tpd_trait <- results$TPDs

# saveRDS(pca_trait, "output/PCA_fish.rds")
# saveRDS(tpd_trait, "output/TPDs_fish.rds")

pca_trait <- readRDS("output/PCA_fish.rds")
tpd_trait <- readRDS("output/TPDs_fish.rds")

################################################################################
# OTHER FISH INFORMATION
################################################################################

fishOtherInfo <- fishData[, (max(columnsTraits) + 1):ncol(fishData)]

################################################################################
# ADD HUMAN USES TO PCA OBJECT
################################################################################

df_scraping <- readRDS("output/fish_human_uses_binary_FB.rds") %>%
  data.table::as.data.table()
df_uni <- data.table::fread("dataPrepared/Fish/uni.csv")

data.table::setnames(
  df_scraping,
  c("species_name", "aquarium", "fisheries", "bait", "game_fish", "aquaculture"),
  c("Species", "Aquarium", "Fisheries", "Bait", "Game_fish", "Aquaculture")
)
data.table::setnames(df_uni, "Game fish", "Game_fish", skip_absent = TRUE)

binary_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game_fish", "Bait")

merged_df <- merge(
  df_uni,
  df_scraping[, c("Species", binary_cols), with = FALSE],
  by = "Species",
  suffixes = c("", ".scraping"),
  all.x = TRUE
)

for (col in binary_cols) {
  col_scraping <- paste0(col, ".scraping")
  if (col_scraping %in% colnames(merged_df)) {
    merged_df[[col]] <- pmax(
      as.numeric(merged_df[[col]]),
      as.numeric(merged_df[[col_scraping]]),
      na.rm = TRUE
    )
    merged_df[[col_scraping]] <- NULL
  }
}

merged_df[, All_uses := as.integer(Fisheries + Aquaculture + Aquarium + Game_fish + Bait > 0)]

merged_df_clean <- merged_df %>%
  dplyr::select(Species, dplyr::all_of(binary_cols), All_uses) %>%
  dplyr::rename("Game fish" = Game_fish, "All uses" = All_uses)

usage_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "Bait", "All uses")

uses_df <- as.data.frame(merged_df_clean[, c("Species", usage_cols), with = FALSE])
uses_df <- uses_df[!is.na(uses_df$Species) & uses_df$Species != "NA", ]
rownames(uses_df) <- uses_df$Species
uses_df$Species <- NULL

# Manual fix for a mismatched name
uses_df["Centromochlus musaica", ] <- list(
  Fisheries = 0,
  Aquaculture = 0,
  Aquarium = 1,
  `Game fish` = 0,
  Bait = 0,
  `All uses` = 1
)

species_pca   <- rownames(pca_trait$traits_scores)
species_uses  <- rownames(uses_df)
species_ref   <- intersect(species_pca, species_uses)

pca_trait$uses <- uses_df[species_ref, usage_cols, drop = FALSE]
use <- pca_trait$uses

# saveRDS(pca_trait, "output/pca_trait.rds")
pca_trait <- readRDS("output/pca_trait.rds")

################################################################################
# COMMUNITY MATRIX (USES x SPECIES)
################################################################################

species_scores <- as.data.frame(pca_trait$traits_scores[, 1:4])
species_uses   <- as.data.frame(pca_trait$uses)

species_uses$rownames   <- rownames(species_uses)
species_scores$rownames <- rownames(species_scores)

species_scores_uses <- merge(species_uses, species_scores, by = "rownames")
rownames(species_scores_uses) <- species_scores_uses$rownames
species_scores_uses$rownames  <- NULL
species_uses$rownames         <- NULL

MatriceFish <- t(as.matrix(species_uses))
column_names <- rownames(species_uses)

MatriceFish_1 <- matrix(1, nrow = 1, ncol = length(column_names))
colnames(MatriceFish_1) <- column_names
rownames(MatriceFish_1) <- "all"

MatriceFish <- rbind(MatriceFish, MatriceFish_1)

# write.csv(MatriceFish, "output/MatriceFish.csv", row.names = TRUE)