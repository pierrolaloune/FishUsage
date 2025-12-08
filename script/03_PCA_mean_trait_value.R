# ============================================================
# Script: 03_PCA_mean_trait_value.R
# Purpose: Null model for PCA mean trait values per human use,
#          SES computation and visualization (PCA distributions
#          and PCA loadings)
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

tpd_trait   <- readRDS("output/TPDs_fish.rds")
pca_trait   <- readRDS("output/PCA_fish.rds")
MatriceFish <- read.csv("output/MatriceFish.csv")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

################################################################################
# NULL MODEL OF PCA MEAN TRAIT VALUE
################################################################################

resultats_null <- generate_null_means( # long recommended skip
  pca_trait     = pca_trait,
  MatriceFish   = MatriceFish,
  nb_simulations = 999
)

resultats_null$all <- NULL

# saveRDS(resultats_null, "output/PCA_mean_trait_values_results.rds")
PCA_mean_trait_values_results <- readRDS("output/PCA_mean_trait_values_results.rds")

################################################################################
# SES COMPUTATION
################################################################################

PCA_mean_trait_values_SES <- get_SES_from_PCA_results( # Table S3
  results_list = PCA_mean_trait_values_results
)

# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")

################################################################################
# BUILD LONG FORMAT FOR PLOTTING
################################################################################

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

################################################################################
# PLOT: SIMULATED PCA DISTRIBUTIONS
################################################################################

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
    strip.text      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5)
  )

print(p)
# ggsave("plot/PCA_mean_trait_values_plot.png", p, width = 6, height = 10, dpi = 300)

################################################################################
# PCA LOADING STRUCTURE
################################################################################

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
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(plot_trait_contrib)
# ggsave("plot/plot_trait_contrib.png", plot_trait_contrib, width = 5, height = 5, dpi = 300)

################################################################################
# PCA CORRELATION CIRCLE (PC1–PC2)
################################################################################

sdev <- pca_trait$pca_object$sdev
var_explained <- sdev^2 / sum(sdev^2)

percent_PC1 <- round(var_explained[1] * 100, 1)
percent_PC2 <- round(var_explained[2] * 100, 1)

x_lab <- paste0("PC1 (", percent_PC1, "%)")
y_lab <- paste0("PC2 (", percent_PC2, "%)")

loadings_df <- as.data.frame(
  pca_trait$pca_object$loadings[, c("Comp.1", "Comp.2")]
)
loadings_df$Trait <- rownames(loadings_df)

circle_df <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 100)),
  y = sin(seq(0, 2 * pi, length.out = 100))
)

cercle_PCA <- ggplot2::ggplot() +
  ggplot2::geom_path(data = circle_df, aes(x = x, y = y), color = "grey70") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  ggplot2::geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "darkred",
    linewidth = 0.3
  ) +
  ggrepel::geom_text_repel(
    data = loadings_df,
    aes(x = Comp.1, y = Comp.2, label = Trait),
    color = "black", size = 3.5, max.overlaps = Inf
  ) +
  ggplot2::coord_fixed() +
  ggplot2::xlim(c(-1.1, 1.1)) +
  ggplot2::ylim(c(-1.1, 1.1)) +
  ggplot2::labs(
    title = "Correlation circle (PCA)",
    x = x_lab,
    y = y_lab
  ) +
  ggplot2::theme_minimal()

print(cercle_PCA)
# ggsave("plot/cercle_PCA.png", cercle_PCA, width = 5, height = 5, dpi = 300)
