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
# simulated PCA distributions with observed values overlaid, (ii) a heatmap of
# signed trait loadings across PCA axes, and (iii) a PCA correlation circle for
# PC1–PC2 based on loadings and explained variance.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait     <- readRDS("output/TPDs_fish.rds")
pca_trait     <- readRDS("output/PCA_fish.rds")
pca_trait_use <- readRDS("output/pca_trait.rds")

MatriceFish <- read.csv("output/MatriceFish.csv")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# ------------------------------------------------------------------------------
# Null model of PCA mean trait value
# ------------------------------------------------------------------------------

resultats_null <- generate_null_means( # long recommended skip
  pca_trait      = pca_trait,
  MatriceFish    = MatriceFish,
  nb_simulations = 999
)

resultats_null$all <- NULL

# saveRDS(resultats_null, "output/PCA_mean_trait_values_results.rds")
PCA_mean_trait_values_results <- readRDS("output/PCA_mean_trait_values_results.rds")

# ------------------------------------------------------------------------------
# SES computation
# ------------------------------------------------------------------------------

PCA_mean_trait_values_SES <- get_SES_from_PCA_results( # Table S3
  results_list = PCA_mean_trait_values_results
)

# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")

# ------------------------------------------------------------------------------
# Build long format for plotting
# ------------------------------------------------------------------------------

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
# ggsave("plot/PCA_mean_trait_values_plot.png", p, width = 6, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# PCA loading structure
# ------------------------------------------------------------------------------

loadings_mat <- as.matrix(
  pca_trait$pca_object$loadings[, paste0("Comp.", 1:4)]
)

contrib_df <- as.data.frame(loadings_mat)
contrib_df$Trait <- rownames(contrib_df)

contrib_long <- contrib_df %>%
  tidyr::pivot_longer(
    cols = starts_with("Comp."),
    names_to = "Component",
    values_to = "Loading"
  ) %>%
  dplyr::mutate(
    Sign         = ifelse(Loading >= 0, "Positive", "Negative"),
    Contribution = 100 * (Loading^2)
  )

plot_trait_contrib <- ggplot2::ggplot(
  contrib_long,
  aes(x = Component, y = reorder(Trait, abs(Loading)), fill = Loading)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    aes(label = sprintf("%.2f", Loading)),
    size = 3.2
  ) +
  ggplot2::scale_fill_gradient2(
    low = "firebrick", mid = "white", high = "steelblue",
    midpoint = 0, name = "Loading\n(coef)"
  ) +
  ggplot2::labs(
    title = "Signed loadings of traits on PCA axes",
    x = "Principal Component",
    y = "Trait"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid  = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

print(plot_trait_contrib)
# ggsave("plot/plot_trait_contrib.png", plot_trait_contrib, width = 5, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(PCA_mean_trait_values_results, "output/PCA_mean_trait_values_results.rds")
# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")
# ggsave("plot/PCA_mean_trait_values_plot.png", p, width = 6, height = 10, dpi = 300)
# ggsave("plot/plot_trait_contrib.png", plot_trait_contrib, width = 5, height = 5, dpi = 300)
