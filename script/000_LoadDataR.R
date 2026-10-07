# --- 000_LoadDataR ------------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

trait      <- readRDS("dataPrepared/Fish/FishMORPH_Traits.rds")
phylogeny  <- readRDS("dataPrepared/Fish/FishMORPH_Phylogeny.rds")
traitNames <- colnames(trait)[-c(1:6)]

# --- FishBase length/weight [LONG] --------------------------------------------

# list_sp <- gsub("\\.", " ", as.character(trait$Genus.species))
#
# lgtwgt <- as.data.frame(rfishbase::length_weight(list_sp))
# lgtwgt <- lgtwgt[!is.na(lgtwgt$a) & !is.na(lgtwgt$b), ]
#
# ab <- lgtwgt %>%
#   dplyr::group_by(Species) %>%
#   dplyr::slice_max(order_by = CoeffDetermination, with_ties = FALSE, na_rm = TRUE) %>%
#   dplyr::select(Species, a, b) %>%
#   dplyr::distinct()
#
# speciesInfo <- rfishbase::species(list_sp) %>% data.table::data.table()
# speciesInfoSub <- unique(speciesInfo[, .(Species, Fresh, LongevityWild, Length, Weight, LTypeMaxM)])
# speciesInfoSub$Species <- gsub(" ", ".", speciesInfoSub$Species)
#
# speciesInfoSub <- speciesInfoSub %>%
#   dplyr::mutate(
#     Length = ifelse(LTypeMaxM != "SL", NA, Length),
#     Weight2 = ifelse(
#       Species %in% gsub(" ", ".", rownames(ab)),
#       ab[gsub("\\.", " ", Species), "a"] * Length^ab[gsub("\\.", " ", Species), "b"],
#       NA
#     )
#   ) %>%
#   dplyr::mutate(Weight = ifelse(is.na(Weight) & !is.na(Weight2), Weight2, Weight))
#
# fishTraits <- merge(
#   trait[, 6:15],
#   speciesInfoSub[, .(Species, Length, Weight)],
#   by.x = "Genus.species", by.y = "Species", all.x = TRUE
# )
#
# fishTraits <- fishTraits %>%
#   dplyr::mutate(dplyr::across(-Genus.species, ~ log10(.x + 1))) %>%
#   dplyr::rename(species = Genus.species) %>%
#   dplyr::mutate(species = gsub("\\.", "_", species))
#
# spToKeep <- fishTraits %>%
#   dplyr::select(species) %>%
#   dplyr::mutate(species = gsub("_", " ", species)) %>%
#   dplyr::pull() %>%
#   rfishbase::species() %>%
#   data.table::as.data.table() %>%
#   dplyr::filter(Fresh == 1) %>%
#   dplyr::select(Species) %>%
#   dplyr::mutate(Species = gsub(" ", "_", Species)) %>%
#   dplyr::pull()
#
# fishTraits <- fishTraits %>% dplyr::filter(species %in% spToKeep)
#
# write.table(fishTraits, "dataPrepared/Fish/fishTraitsMissing.txt")

fishTraits <- read.table("dataPrepared/Fish/fishTraitsMissing.txt")

# --- Phylogenetic PCoA [LONG] -------------------------------------------------

# phylogeny$tip.label <- gsub("\\.", "_", phylogeny$tip.label)
# phylogeny <- ape::drop.tip(phylogeny, setdiff(phylogeny$tip.label, fishTraits$species))
# phylogenyTraits <- phytools::force.ultrametric(phylogeny)
# phylDiss <- sqrt(cophenetic(phylogenyTraits))
#
# pcoaPhyl <- cmdscale(phylDiss, k = 10)
# write.table(pcoaPhyl, "dataPrepared/Fish/pcoaPhylogenyFish.txt")

pcoaPhyl <- read.table("dataPrepared/Fish/pcoaPhylogenyFish.txt", header = TRUE, stringsAsFactors = FALSE)

rownames(pcoaPhyl) <- fishTraits$species
colnames(pcoaPhyl) <- paste0("Eigen.", 1:10)

# --- Taxonomy [LONG] ----------------------------------------------------------

# list_sp_raw <- stringr::str_squish(gsub("_", " ", fishTraits$species))
#
# gna_one_safe <- function(x) {
#   tryCatch(
#     taxize::gna_verifier(
#       names = x, data_sources = 11, all_matches = FALSE,
#       capitalize = TRUE, species_group = TRUE, output_type = "table"
#     ),
#     error = function(e) {
#       tibble::tibble(submittedName = x, matchedName = NA_character_,
#                      matchType = "Error", dataSourceId = NA_real_)
#     }
#   )
# }
#
# progressr::handlers(global = TRUE)
# progressr::handlers("txtprogressbar")
#
# verified_names <- progressr::with_progress({
#   p <- progressr::progressor(along = list_sp_raw)
#   purrr::map_dfr(list_sp_raw, function(x) {
#     p(message = x)
#     gna_one_safe(x)
#   })
# })
#
# new_names <- gsub(" ", "_", verified_names$matchedCanonicalSimple)
#
# rownames(pcoaPhyl) <- new_names
# fishTraits$species <- new_names
# traitsAndPCOA      <- cbind(fishTraits, pcoaPhyl)
#
# write.table(traitsAndPCOA, "dataPrepared/Fish/traitsWithPCOA.txt")

traitsAndPCOA <- read.table("dataPrepared/Fish/traitsWithPCOA.txt", header = TRUE, stringsAsFactors = FALSE)

# --- IUCN [LONG] --------------------------------------------------------------

# species_list        <- unique(stringr::str_trim(traitsAndPCOA$species))
# species_list_spaces <- gsub("_", " ", species_list)
#
# iucn_clean <- read.csv("dataOriginal/assessments.csv", header = TRUE, stringsAsFactors = FALSE) %>%
#   dplyr::mutate(scientificName = stringr::str_trim(scientificName))
#
# species_to_update <- tibble::tibble(species = species_list_spaces) %>%
#   dplyr::anti_join(iucn_clean, by = c("species" = "scientificName"))
#
# synonyms_mapping <- rfishbase::synonyms(
#   species_list = species_to_update$species,
#   server = "fishbase", version = "latest", fields = NULL
# ) %>%
#   dplyr::filter(!Status %in% c("misapplied name", "ambiguous synonym", "provisionally accepted name")) %>%
#   dplyr::select(Species, synonym) %>%
#   dplyr::distinct()
#
# iucn_clean <- iucn_clean %>%
#   dplyr::mutate(
#     scientificName = dplyr::if_else(
#       scientificName %in% synonyms_mapping$Species,
#       synonyms_mapping$synonym[match(scientificName, synonyms_mapping$Species)],
#       scientificName
#     )
#   )
#
# species_to_check <- read.csv("dataOriginal/species_to_update_900_done.csv",
#                              sep = ";", header = TRUE, stringsAsFactors = FALSE)
# colnames(species_to_check) <- c("scientificName", "redlistCategory")
#
# acronyms <- c(
#   "Critically Endangered" = "CR", "Endangered" = "EN", "Vulnerable" = "VU",
#   "Near Threatened" = "NT", "Least Concern" = "LC", "Data Deficient" = "DD",
#   "Extinct" = "EX", "Extinct in the Wild" = "EW", "Not Evaluated" = "NE"
# )
#
# iucn_clean <- iucn_clean[, c("scientificName", "redlistCategory")] %>%
#   dplyr::mutate(redlistCategory = dplyr::recode(redlistCategory, !!!acronyms)) %>%
#   dplyr::bind_rows(species_to_check) %>%
#   dplyr::filter(scientificName %in% species_list_spaces) %>%
#   dplyr::group_by(scientificName) %>%
#   dplyr::slice(1) %>%
#   dplyr::ungroup()
#
# traitsAndPCOA$species <- gsub("_", " ", traitsAndPCOA$species)
# traitsAndPCOA$IUCN <- iucn_clean$redlistCategory[match(traitsAndPCOA$species, iucn_clean$scientificName)]
#
# write.table(traitsAndPCOA, "dataPrepared/Fish/traitsWithPCOAIUCN.txt")

fishData <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt", header = TRUE, stringsAsFactors = FALSE)

# --- missForest imputation [LONG] ---------------------------------------------

columnsTraits     <- 2:(which(colnames(fishData) == "Eigen.1") - 1)
columnsImputation <- 2:(which(colnames(fishData) == "IUCN") - 1)

# set.seed(123)
# imputed_forest <- missForest::missForest(xmis = fishData[, columnsImputation])
# print(imputed_forest$OOBerror)
#
# fishData_imputed_forest <- fishData
# fishData_imputed_forest[colnames(fishData)[columnsTraits]] <- imputed_forest$ximp[colnames(fishData)[columnsTraits]]
#
# write.table(fishData_imputed_forest, "dataPrepared/Fish/fishData_imputed_forest.txt", row.names = FALSE)

fishData_imputed_forest <- read.table("dataPrepared/Fish/fishData_imputed_forest.txt", header = TRUE)

# --- Trait selection ----------------------------------------------------------

selectedTraits <- c("EdHd", "EhBd", "JlHd", "MoBd", "BlBd", "HdBd",
                    "PFiBd", "PFlBl", "CFdCPd", "Length", "Weight")
newTraitNames  <- c("es", "ep", "ms", "mp", "elo", "wid",
                    "pp", "ps", "cs", "svl", "bm")

fishTraitsMissing <- fishData[, selectedTraits]
fishTraitsImputed <- fishData_imputed_forest[, selectedTraits]
colnames(fishTraitsMissing) <- colnames(fishTraitsImputed) <- newTraitNames
rownames(fishTraitsMissing) <- rownames(fishTraitsImputed) <- fishData$species

fishTraitsMissing <- data.frame(fishTraitsMissing, IUCN = fishData$IUCN)
fishTraitsImputed <- data.frame(fishTraitsImputed, IUCN = fishData$IUCN)

# write.table(fishTraitsMissing, "dataPrepared/Fish/TraitFishMissing.txt")
# write.table(fishTraitsImputed, "dataPrepared/Fish/TraitFishImputed.txt")

# --- PCA + TPD [LONG] ---------------------------------------------------------

# results <- computePCAandTPDs(fishTraitsImputed[, colnames(fishTraitsImputed) != "IUCN"])
# saveRDS(results, "output/All_fish.rds")
# saveRDS(results$PCA,  "output/PCA_fish.rds")
# saveRDS(results$TPDs, "output/TPDs_fish.rds")

pca_trait <- readRDS("output/PCA_fish.rds")
tpd_trait <- readRDS("output/TPDs_fish.rds")

# --- Human uses ---------------------------------------------------------------

df_scraping <- data.table::as.data.table(readRDS("output/fish_human_uses_binary_FB.rds"))
df_uni      <- data.table::fread("dataPrepared/Fish/uni.csv")

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
  by = "Species", suffixes = c("", ".scraping"), all.x = TRUE
)

for (col in binary_cols) {
  col_scraping <- paste0(col, ".scraping")
  merged_df[[col]] <- pmax(as.numeric(merged_df[[col]]), as.numeric(merged_df[[col_scraping]]), na.rm = TRUE)
  merged_df[[col_scraping]] <- NULL
}

merged_df[, All_uses := as.integer(Fisheries + Aquaculture + Aquarium + Game_fish + Bait > 0)]

usage_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "All uses")

uses_df <- merged_df %>%
  dplyr::select(Species, dplyr::all_of(binary_cols), All_uses) %>%
  dplyr::rename("Game fish" = Game_fish, "All uses" = All_uses) %>%
  as.data.frame()
uses_df <- uses_df[!is.na(uses_df$Species) & uses_df$Species != "NA", c("Species", usage_cols)]
rownames(uses_df) <- uses_df$Species
uses_df$Species <- NULL

uses_df["Centromochlus musaica", ] <- list(0, 0, 1, 0, 1)

species_ref    <- intersect(rownames(pca_trait$traits_scores), rownames(uses_df))
pca_trait$uses <- uses_df[species_ref, usage_cols, drop = FALSE]

# saveRDS(pca_trait, "output/pca_trait.rds")
pca_trait <- readRDS("output/pca_trait.rds")

# --- Community matrix ---------------------------------------------------------

MatriceFish <- rbind(
  t(as.matrix(pca_trait$uses)),
  all = rep(1, nrow(pca_trait$uses))
)

# write.csv(MatriceFish, "output/MatriceFish.csv", row.names = TRUE)

# --- Session ------------------------------------------------------------------

sessionInfo()
