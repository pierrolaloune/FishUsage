# ============================================================
# Script: 000_ScrappingData.R
# Purpose: Scrape human use information from FishBase
#          and generate binary use variables per species
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

uni <- data.table::fread("dataPrepared/Fish/uni.csv")
data.table::setnames(uni, old = colnames(uni)[8], new = "Non_uses")

species_names <- uni$Species

################################################################################
# SCRAPING FROM FISHBASE (PARALLEL + PROGRESS)
################################################################################

# recommended not run it, it's long
# Build FishBase URLs for each species
# species_urls <- purrr::map_chr(species_names, make_fishbase_url)
# 
# future::plan(future::multisession)
# progressr::handlers(global = TRUE)
# 
# results_raw <- progressr::with_progress({
#   p <- progressr::progressor(along = species_urls)
#   furrr::future_map_dfr(species_urls, ~ {
#     p()
#     extract_human_uses(.x)
#   })
# })
# 
# future::plan(future::sequential)
# 
# results_raw <- results_raw %>%
#   dplyr::mutate(species_name = species_names)

# saveRDS(results_raw, file = "output/results_raw.rds")

results_raw <- readRDS("output/results_raw.rds")

################################################################################
# CLASSIFICATION OF USE TYPES
################################################################################

results_classified <- results_raw %>%
  dplyr::bind_cols(
    purrr::map_dfr(results_raw$human_uses, classify_uses_precise)
  )

results_classified_selected <- results_classified %>%
  dplyr::select(species_name, aquarium, fisheries, bait, game_fish, aquaculture) %>%
  data.table::as.data.table()

#saveRDS(results_classified_selected, "output/fish_human_uses_classified_parallel_progress.rds")
results_classified_selected <- readRDS("output/fish_human_uses_classified_parallel_progress.rds")

# Example check:
# results_classified %>% dplyr::filter(species_name == "Abramis brama") %>% as.data.frame()

################################################################################
# FINAL BINARY VARIABLES
################################################################################

results_binary <- results_classified_selected %>%
  dplyr::mutate(dplyr::across(
    .cols = -species_name,
    .fns = ~ dplyr::case_when(
      .x == "none" ~ 0L,
      .x == "rare" ~ 0L,
      is.na(.x)    ~ 0L,
      TRUE         ~ 1L
    )
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(all_use = if_else(sum(dplyr::c_across(-species_name)) > 0, 1L, 0L)) %>%
  dplyr::ungroup()

# saveRDS(results_binary, "output/fish_human_uses_binary_FB.rds")
