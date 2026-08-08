# ------------------------------------------------------------------------------
# Script : web scrapping percentage
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script measures how much the FishBase web scraping actually added to the
# use data, on top of what the rfishbase tables already provided.
#
# For each use, and then for all uses together, the species are split into three
# groups: known from rfishbase only, known from scraping only, and known from
# both. The percentages printed at the end are the ones quoted in the paper.
#
# The script ends with a second, finer pass that works on merged_df, the merged
# table built in 000_LoadDataR.R, and reports the same split per use.
#
# NOTE: this script does not create merged_df. Run 000_LoadDataR.R first, or the
# second half will fail with "object 'merged_df' not found".

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

df_scraping <- readRDS("output/fish_human_uses_binary_FB.rds")
df_uni <- read.csv("dataPrepared/Fish/uni.csv")
pca_trait <- readRDS("output/pca_trait.rds")
species_ref <- rownames(pca_trait$traits_scaled)

# ---- Align the column names of the two sources so they can be stacked ----
colnames(df_scraping) <- c("Species", "Aquarium", "Fisheries", "Bait",
                           "Game_fish", "Aquaculture", "all_use")
colnames(df_uni)[colnames(df_uni) == "Game.fish"] <- "Game_fish"

# ---- Add one species missing from both tables (Centromochlus musaica) ----
new_row_scrap <- data.table(Species = "Centromochlus musaica", Aquarium = 1,
                            Fisheries = 0, Bait = 0, Game_fish = 0,
                            Aquaculture = 0, all_use = 1)

new_row_uni <- data.table(Species = "Centromochlus musaica", Fisheries = 0,
                          Aquaculture = 0, Aquarium = 1, Game_fish = 0,
                          Bait = 0, Ui = NA, `Non uses` = 0)

df_scraping <- rbindlist(list(df_scraping, new_row_scrap), use.names=TRUE, fill=TRUE)
df_uni <- rbindlist(list(df_uni, new_row_uni), use.names=TRUE, fill=TRUE)

# ------------------------------------------------------------------------------
# Filter and align on the reference species list
# ------------------------------------------------------------------------------

# Only the species kept in the PCA are counted, so that both sources describe
# exactly the same pool.
df_uni_clean   <- df_uni[Species %in% species_ref]
df_scrap_clean <- df_scraping[Species %in% species_ref]

# Sort by species, so that row i of one table matches row i of the other.
setkey(df_uni_clean, Species)
setkey(df_scrap_clean, Species)

# ------------------------------------------------------------------------------
# Counts and percentages per use
# ------------------------------------------------------------------------------

target_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game_fish")
results <- list()

for (usage in target_cols) {
  # Logical vector per source: is this species recorded for this use?
  val_uni   <- df_uni_clean[[usage]] == 1
  val_scrap <- df_scrap_clean[[usage]] == 1

  # Exclusive segments
  count_uni_only   <- sum(val_uni & !val_scrap, na.rm = TRUE)
  count_scrap_only <- sum(!val_uni & val_scrap, na.rm = TRUE) # net gain of the scraping
  count_both       <- sum(val_uni & val_scrap, na.rm = TRUE)

  total_unique_species <- sum(val_uni | val_scrap, na.rm = TRUE)

  results[[usage]] <- data.table(
    Usage              = usage,
    `rFishBase Only`   = count_uni_only,
    `WebScraping Only` = count_scrap_only,
    `Shared (Merge)`   = count_both,
    `Total Unique`     = total_unique_species,
    `% Scraping Gain`  = round((count_scrap_only / total_unique_species) * 100, 2)
  )
}

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------

df_final_report <- rbindlist(results)

print(df_final_report)

# ------------------------------------------------------------------------------
# Global figures (all species, all uses pooled)
# ------------------------------------------------------------------------------

# A species counts as documented as soon as it has at least one of the four uses.
global_uni   <- rowSums(df_uni_clean[, ..target_cols] == 1, na.rm = TRUE) > 0
global_scrap <- rowSums(df_scrap_clean[, ..target_cols] == 1, na.rm = TRUE) > 0

# Overall counts (full Venn diagram)
n_total_species_with_info <- sum(global_uni | global_scrap)
n_glob_uni_only   <- sum(global_uni & !global_scrap)
n_glob_scrap_only <- sum(!global_uni & global_scrap)
n_glob_both       <- sum(global_uni & global_scrap)

# Percentages quoted in the manuscript
p_glob_uni_only   <- round((n_glob_uni_only / n_total_species_with_info) * 100, 1)
p_glob_scrap_only <- round((n_glob_scrap_only / n_total_species_with_info) * 100, 1)
p_glob_both       <- round((n_glob_both / n_total_species_with_info) * 100, 1)

cat("\n--- Values for the manuscript ---\n")
cat("Total species with info:", n_total_species_with_info, "\n")
cat("Exclusively pre-assembled:", p_glob_uni_only, "%\n")
cat("Exclusively scraping:", p_glob_scrap_only, "%\n")
cat("Both sources:", p_glob_both, "%\n")


# ---- Same split, recomputed from the summary table above ----
total_rfishbase <- sum(df_final_report$`rFishBase Only`)
total_webscraping <- sum(df_final_report$`WebScraping Only`)
total_shared <- sum(df_final_report$`Shared (Merge)`)
total_unique <- sum(df_final_report$`Total Unique`)

cat("=== GLOBAL SUMMARY OF THE HUMAN USE DATA ===\n\n")
cat(sprintf("Total rFishBase Only:     %d (%.1f%%)\n",
            total_rfishbase, 100 * total_rfishbase / total_unique))
cat(sprintf("Total WebScraping Only:   %d (%.1f%%)\n",
            total_webscraping, 100 * total_webscraping / total_unique))
cat(sprintf("Total Shared (Merge):     %d (%.1f%%)\n",
            total_shared, 100 * total_shared / total_unique))
cat(sprintf("Total Unique:             %d (100.0%%)\n\n", total_unique))

cat(sprintf("Global scraping gain: %.1f%%\n",
            100 * total_webscraping / (total_rfishbase + total_unique)))

rfishbase_pct <- 100 * total_rfishbase / total_unique
webscraping_pct <- 100 * total_webscraping / total_unique
shared_pct <- 100 * total_shared / total_unique

cat(sprintf("\n-> from rFishBase only:   %.1f%%\n", rfishbase_pct))
cat(sprintf("-> from web scraping only: %.1f%%\n", webscraping_pct))
cat(sprintf("-> from both sources:      %.1f%%\n", shared_pct))



# ------------------------------------------------------------------------------
# Second pass: same split computed on merged_df (from 000_LoadDataR.R)
# ------------------------------------------------------------------------------

df_scraping <- readRDS("output/fish_human_uses_binary_FB.rds") %>%
  data.table::as.data.table()
df_uni <- data.table::fread("dataPrepared/Fish/uni.csv") # from rfishbase
pca_trait  <- readRDS("output/pca_trait.rds")
species_ref <- rownames(pca_trait$traits_scaled)
ref_species_list <- as.character(species_ref)
df_scraping <- df_scraping[species_name  %in% ref_species_list]
merged_df <- merged_df[Species %in% ref_species_list]

library(data.table)

# ---- Filter both tables on the reference species list ----
ref_species_list <- as.character(species_ref)

df_scraping_filtered <- df_scraping[species_name %in% ref_species_list]
merged_df_filtered <- merged_df[Species %in% ref_species_list]

# ---- Keep only the species that have at least one documented use ----
usages <- c("Fisheries", "Aquaculture", "Aquarium", "Game_fish")
usages_scrap <- paste0(usages, ".scraping")

# A row is informative if it holds at least one 1, in either source.
has_any_use <- rowSums(merged_df_filtered[, c(..usages, ..usages_scrap), with = FALSE] == 1, na.rm = TRUE) > 0

df_analysis <- merged_df_filtered[has_any_use]
n_species_analysis <- nrow(df_analysis)

# ---- Recap table, one row per use ----
source_summary <- data.table(
  Usage = usages,
  Total_Final_1 = 0,
  Rfishbase_only_pct = 0,
  Scraping_only_pct = 0,
  Both_pct = 0
)

for (i in seq_along(usages)) {
  col <- usages[i]
  col_scrap <- paste0(col, ".scraping")

  if (col_scrap %in% colnames(df_analysis)) {

    # Extract both vectors, guarding against NA and type mismatches
    r_val <- as.numeric(df_analysis[[col]] == 1 & !is.na(df_analysis[[col]]))
    s_val <- as.numeric(df_analysis[[col_scrap]] == 1 & !is.na(df_analysis[[col_scrap]]))

    # Intersections for this specific use
    total_1 <- sum(r_val == 1 | s_val == 1)

    if (total_1 > 0) {
      r_only <- sum(r_val == 1 & s_val == 0)
      s_only <- sum(r_val == 0 & s_val == 1)
      both   <- sum(r_val == 1 & s_val == 1)

      source_summary$Total_Final_1[i]      <- total_1
      source_summary$Rfishbase_only_pct[i] <- round((r_only / total_1) * 100, 1)
      source_summary$Scraping_only_pct[i]  <- round((s_only / total_1) * 100, 1)
      source_summary$Both_pct[i]           <- round((both / total_1) * 100, 1)
    }
  }
}

# ------------------------------------------------------------------------------
# Sentence used in the manuscript
# ------------------------------------------------------------------------------

total_all_records <- sum(source_summary$Total_Final_1)

# Global averages, weighted by the number of species per use
global_rf_only <- round(sum(source_summary$Total_Final_1 * source_summary$Rfishbase_only_pct / 100) / total_all_records * 100, 1)
global_sc_only <- round(sum(source_summary$Total_Final_1 * source_summary$Scraping_only_pct / 100) / total_all_records * 100, 1)
global_both    <- round(sum(source_summary$Total_Final_1 * source_summary$Both_pct / 100) / total_all_records * 100, 1)

print(source_summary)

cat(sprintf(
  "\nAmong the %s species for which use information was found, %s%% were recorded exclusively from HTML FishBase pages, %s%% exclusively by scraping, and %s%% by both sources.",
  format(n_species_analysis, big.mark=","), global_rf_only, global_sc_only, global_both
))
