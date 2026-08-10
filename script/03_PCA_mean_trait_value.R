# ------------------------------------------------------------------------------
# Script : 03_PCA_mean_trait_value
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads fish TPDs, PCA outputs, and the community-by-use matrix, then
# runs a null model to simulate mean PCA scores per human-use community. It
# computes SES for observed mean PCA values against simulated distributions and
# prepares long-format data for visualization. It produces (i) histograms of
# simulated PCA distributions with observed values overlaid, and (ii) a heatmap
# of signed trait loadings across PCA axes.
#
# Output: Table S3 (SES of the mean PCA position per use).
#
# The null model is commented out; its result is stored in output/ and reloaded
# immediately afterwards, so the script runs end to end as it is.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait     <- readRDS("output/TPDs_fish.rds")
pca_trait     <- readRDS("output/PCA_fish.rds")
pca_trait_use <- readRDS("output/pca_trait.rds")

MatriceFish <- read.csv("output/MatriceFish.csv")

# ---- Restore the species names and the row names lost by read.csv ----
colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# ------------------------------------------------------------------------------
# Null model of PCA mean trait value  [LONG]
# ------------------------------------------------------------------------------

# LONG: 999 randomizations of the use assignments, per human use.

# resultats_null <- generate_null_means(
#   pca_trait      = pca_trait,
#   MatriceFish    = MatriceFish,
#   nb_simulations = 999
# )
#
# # Drop the global "all" category, which has no null distribution of its own.
# resultats_null$all <- NULL
#
# saveRDS(resultats_null, "output/PCA_mean_trait_values_results.rds")

# ---- Recommended: load the saved result ----
PCA_mean_trait_values_results <- readRDS("output/PCA_mean_trait_values_results.rds")
for (usage in names(PCA_mean_trait_values_results)) {
  PCA_mean_trait_values_results[[usage]]$observed["Comp.1"] <-
    -PCA_mean_trait_values_results[[usage]]$observed["Comp.1"]
  
  PCA_mean_trait_values_results[[usage]]$simulated[, "Comp.1"] <-
    -PCA_mean_trait_values_results[[usage]]$simulated[, "Comp.1"]
}
# ------------------------------------------------------------------------------
# SES computation
# ------------------------------------------------------------------------------

# ---- Table S3 ----
PCA_mean_trait_values_SES <- get_SES_from_PCA_results(
  results_list = PCA_mean_trait_values_results
)

# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")

# ------------------------------------------------------------------------------
# Build long format for plotting
# ------------------------------------------------------------------------------

# One row per simulated value, with the matching observed value repeated so that
# ggplot can draw both on the same facet.

df_plot <- purrr::map_dfr(names(PCA_mean_trait_values_results), function(usage) {

  simulated <- as.data.frame(
    PCA_mean_trait_values_results[[usage]]$simulated
  )[, c("Comp.1", "Comp.2")]

  simulated_long <- tidyr::pivot_longer(
    simulated,
    cols      = everything(),
    names_to  = "Component",
    values_to = "Simulated_value"
  )

  observed <- PCA_mean_trait_values_results[[usage]]$observed[c("Comp.1", "Comp.2")]

  simulated_long %>%
    dplyr::mutate(
      Usage          = usage,
      Observed_value = observed[Component]
    )
})

df_plot$Usage <- factor(
  df_plot$Usage,
  levels = c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "Bait", "All uses")
)

# ------------------------------------------------------------------------------
# Plot: simulated PCA distributions
# ------------------------------------------------------------------------------

p <- ggplot2::ggplot(df_plot, aes(x = Simulated_value)) +
  ggplot2::geom_histogram(bins = 50, fill = "#69b3a2", color = "black") +
  ggplot2::geom_vline(
    aes(xintercept = Observed_value),
    color = "red", linetype = "dashed", linewidth = 0.8
  ) +
  ggplot2::facet_grid(Usage ~ Component, scales = "free") +
  ggplot2::labs(
    x = "Simulated PCA score",
    y = "Frequency",
    title = "Simulated PCA distributions vs Observed values",
    subtitle = "PCA Components 1 and 2"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    strip.text    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(PCA_mean_trait_values_results, "output/PCA_mean_trait_values_results.rds")
# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")
# ggsave("plot/PCA_mean_trait_values_plot.png", p, width = 6, height = 10, dpi = 300)
