# ------------------------------------------------------------------------------
# Script : 000_ScrappingData
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads a reference species list, scrapes human-use information from
# FishBase (optionally in parallel with progress reporting), classifies each
# record into standardized use categories, and converts these categories into
# final binary use variables per species. The resulting binary table can be saved
# for downstream analyses.
#
# The scraping step itself is commented out: it takes days, because a random
# 5-8 s pause is applied before each request. Its result is stored in output/
# and reloaded below, so the script runs end to end as it is.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
uni <- data.table::fread("dataPrepared/Fish/uni.csv")
data.table::setnames(uni, old = colnames(uni)[8], new = "Non_uses")

species_names <- uni$Species

# ------------------------------------------------------------------------------
# Scraping from FishBase  [LONG]
# ------------------------------------------------------------------------------

# LONG: one HTTP request per species, with a mandatory pause between requests.
# Result saved as results_raw.rds and reloaded below.

# # ---- Build one FishBase URL per species ----
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
#
# saveRDS(results_raw, file = "output/results_raw.rds")

# ---- Recommended: load the saved result ----
results_raw <- readRDS("output/results_raw.rds")

# ------------------------------------------------------------------------------
# Classification of use types
# ------------------------------------------------------------------------------

# Turn the free-text "Human uses" field into one intensity level per category
# ("none", "rare", "regular", "highly"). See classify_uses_precise() in
# 000_functions.R for the rules.

results_classified <- results_raw %>%
  dplyr::bind_cols(
    purrr::map_dfr(results_raw$human_uses, classify_uses_precise)
  )

results_classified_selected <- results_classified %>%
  dplyr::select(species_name, aquarium, fisheries, bait, game_fish, aquaculture) %>%
  data.table::as.data.table()

# saveRDS(results_classified_selected, "output/fish_human_uses_classified_parallel_progress.rds")
results_classified_selected <- readRDS("output/fish_human_uses_classified_parallel_progress.rds")

# ---- Quick check on a single species ----
results_classified_selected %>% dplyr::filter(species_name == "Abramis brama") %>% as.data.frame()

# ------------------------------------------------------------------------------
# Final binary variables
# ------------------------------------------------------------------------------

# Only "regular" and "highly" count as an actual use: "rare" and "none" are set
# to 0, so a species is flagged as used only when the use is well established.

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

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(results_binary, "output/fish_human_uses_binary_FB.rds")
