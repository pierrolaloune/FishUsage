# --- web scrapping percentage -------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

df_scraping <- readRDS("output/fish_human_uses_binary_FB.rds")
df_uni      <- read.csv("dataPrepared/Fish/uni.csv")
pca_trait   <- readRDS("output/pca_trait.rds")
species_ref <- rownames(pca_trait$traits_scaled)

colnames(df_scraping) <- c("Species", "Aquarium", "Fisheries", "Bait",
                           "Game_fish", "Aquaculture", "all_use")
colnames(df_uni)[colnames(df_uni) == "Game.fish"] <- "Game_fish"

df_scraping <- rbindlist(list(df_scraping, data.table(
  Species = "Centromochlus musaica", Aquarium = 1, Fisheries = 0, Bait = 0,
  Game_fish = 0, Aquaculture = 0, all_use = 1
)), use.names = TRUE, fill = TRUE)
df_uni <- rbindlist(list(df_uni, data.table(
  Species = "Centromochlus musaica", Fisheries = 0, Aquaculture = 0, Aquarium = 1,
  Game_fish = 0, Bait = 0, Ui = NA, `Non uses` = 0
)), use.names = TRUE, fill = TRUE)

df_uni_clean   <- df_uni[Species %in% species_ref]
df_scrap_clean <- df_scraping[Species %in% species_ref]
setkey(df_uni_clean, Species)
setkey(df_scrap_clean, Species)

target_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game_fish")

# --- Sources ------------------------------------------------------------------

source_of_use <- function(usage) {
  val_uni   <- df_uni_clean[[usage]] == 1
  val_scrap <- df_scrap_clean[[usage]] == 1
  dplyr::case_when(
    val_uni & val_scrap  ~ "Both",
    val_uni & !val_scrap ~ "rFishBase only",
    !val_uni & val_scrap ~ "Scraping_FB only"
  )
}

df_final_report <- rbindlist(lapply(target_cols, function(usage) {
  src <- source_of_use(usage)
  data.table(
    Usage              = usage,
    `rFishBase Only`   = sum(src == "rFishBase only", na.rm = TRUE),
    `WebScraping Only` = sum(src == "Scraping_FB only", na.rm = TRUE),
    `Shared (Merge)`   = sum(src == "Both", na.rm = TRUE),
    `Total Unique`     = sum(!is.na(src)),
    `% Scraping Gain`  = round(100 * sum(src == "Scraping_FB only", na.rm = TRUE) / sum(!is.na(src)), 2)
  )
}))

print(df_final_report)

global_uni   <- rowSums(df_uni_clean[, ..target_cols] == 1, na.rm = TRUE) > 0
global_scrap <- rowSums(df_scrap_clean[, ..target_cols] == 1, na.rm = TRUE) > 0

n_total_species_with_info <- sum(global_uni | global_scrap)

cat("\n--- Values for the manuscript ---\n")
cat("Total species with info:", n_total_species_with_info, "\n")
cat("Exclusively pre-assembled:", round(100 * sum(global_uni & !global_scrap) / n_total_species_with_info, 1), "%\n")
cat("Exclusively scraping:", round(100 * sum(!global_uni & global_scrap) / n_total_species_with_info, 1), "%\n")
cat("Both sources:", round(100 * sum(global_uni & global_scrap) / n_total_species_with_info, 1), "%\n")

total_rfishbase   <- sum(df_final_report$`rFishBase Only`)
total_webscraping <- sum(df_final_report$`WebScraping Only`)
total_shared      <- sum(df_final_report$`Shared (Merge)`)
total_unique      <- sum(df_final_report$`Total Unique`)

cat("\n=== GLOBAL SUMMARY OF THE HUMAN USE DATA ===\n")
cat(sprintf("Total rFishBase Only:     %d (%.1f%%)\n", total_rfishbase, 100 * total_rfishbase / total_unique))
cat(sprintf("Total WebScraping Only:   %d (%.1f%%)\n", total_webscraping, 100 * total_webscraping / total_unique))
cat(sprintf("Total Shared (Merge):     %d (%.1f%%)\n", total_shared, 100 * total_shared / total_unique))
cat(sprintf("Total Unique:             %d (100.0%%)\n", total_unique))
cat(sprintf("Global scraping gain: %.1f%%\n", 100 * total_webscraping / (total_rfishbase + total_unique)))

merged_df <- merge(
  as.data.table(data.table::fread("dataPrepared/Fish/uni.csv"))[, Game_fish := `Game fish`],
  as.data.table(readRDS("output/fish_human_uses_binary_FB.rds"))[
    , .(Species = species_name, Fisheries = fisheries, Aquaculture = aquaculture,
        Aquarium = aquarium, Game_fish = game_fish)],
  by = "Species", suffixes = c("", ".scraping"), all.x = TRUE
)[Species %in% species_ref]

usages       <- c("Fisheries", "Aquaculture", "Aquarium", "Game_fish")
usages_scrap <- paste0(usages, ".scraping")

df_analysis <- merged_df[rowSums(merged_df[, c(usages, usages_scrap), with = FALSE] == 1, na.rm = TRUE) > 0]

source_summary <- rbindlist(lapply(usages, function(col) {
  r_val <- df_analysis[[col]] == 1 & !is.na(df_analysis[[col]])
  s_val <- df_analysis[[paste0(col, ".scraping")]] == 1 & !is.na(df_analysis[[paste0(col, ".scraping")]])
  total_1 <- sum(r_val | s_val)
  data.table(
    Usage              = col,
    Total_Final_1      = total_1,
    Rfishbase_only_pct = round(100 * sum(r_val & !s_val) / total_1, 1),
    Scraping_only_pct  = round(100 * sum(!r_val & s_val) / total_1, 1),
    Both_pct           = round(100 * sum(r_val & s_val) / total_1, 1)
  )
}))

print(source_summary)

total_all_records <- sum(source_summary$Total_Final_1)
global_mean <- function(col) round(sum(source_summary$Total_Final_1 * source_summary[[col]] / 100) / total_all_records * 100, 1)

cat(sprintf(
  "\nAmong the %s species for which use information was found, %s%% were recorded exclusively from HTML FishBase pages, %s%% exclusively by scraping, and %s%% by both sources.\n",
  format(nrow(df_analysis), big.mark = ","), global_mean("Rfishbase_only_pct"),
  global_mean("Scraping_only_pct"), global_mean("Both_pct")
))

# --- Session ------------------------------------------------------------------

sessionInfo()
