# ------------------------------------------------------------------------------
# Script : plot_pca_correlation_circle
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script draws the PCA correlation circle of the morphological traits, in
# three independent versions that can be run one after the other:
#   1. Labelled circle - one arrow per trait, coloured by loading magnitude, with
#      the full trait names placed by ggrepel.
#   2. Bare circle - the same arrows in a single colour, without labels, meant to
#      be used as an inset in another figure.
#   3. Heatmap - trait x axis Pearson correlations for PC1 to PC4, with
#      Bonferroni-corrected significance stars.
#
# PC1 is flipped in every version, so the orientation matches the other figures
# of the paper.
#
# Outputs: figures/pca_correlation_circle_highres.png,
#          figures/pca_correlation_circle.png,
#          figures/figS4_clean.png.

# ------------------------------------------------------------------------------
# 1. Labelled correlation circle
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(viridis) # robust, accessible color scales

# ---- Data preparation ----
pca_trait <- readRDS("output/pca_trait.rds")
pca_obj <- pca_trait$pca_object
loadings_raw <- pca_obj$loadings[, 1:2]

# Flip PC1 so that the orientation matches the other figures
loadings_raw[, 1] <- -loadings_raw[, 1]

# Variance explained
sdev          <- pca_obj$sdev
var_explained <- sdev^2 / sum(sdev^2)
pct_pc1       <- round(var_explained[1] * 100, 1)
pct_pc2       <- round(var_explained[2] * 100, 1)

# Full descriptive trait labels
trait_labels <- c(
  es  = "Relative eye size",
  ep  = "Vertical eye position",
  ms  = "Relative maxillary length",
  mp  = "Oral gape position",
  elo = "Body elongation",
  wid = "Body lateral shape",
  pp  = "Pectoral fin vertical position",
  ps  = "Pectoral fin size",
  cs  = "Caudal peduncle throttling",
  svl = "Standard length",
  bm  = "Body mass"
)

loadings_df <- loadings_raw |>
  as.data.frame() |>
  rownames_to_column("trait_abbr") |>
  rename(PC1 = Comp.1, PC2 = Comp.2) |>
  mutate(
    label     = dplyr::recode(trait_abbr, !!!trait_labels),
    magnitude = sqrt(PC1^2 + PC2^2)
  )

# Unit circle data
circle_df <- tibble(
  angle = seq(0, 2 * pi, length.out = 300),
  x     = cos(angle),
  y     = sin(angle)
)

# ---- Visualization ----

# Same table as above, plus a label offset that depends on the quadrant, so that
# the trait names are pushed outwards rather than towards the centre.
loadings_df <- loadings_raw |>
  as.data.frame() |>
  rownames_to_column("trait_abbr") |>
  rename(PC1 = Comp.1, PC2 = Comp.2) |>
  mutate(
    label     = dplyr::recode(trait_abbr, !!!trait_labels),
    magnitude = sqrt(PC1^2 + PC2^2),
    # offset per quadrant
    nudge_x = case_when(
      PC1 > 0 & PC2 > 0 ~  0.08,   # quadrant ++
      PC1 < 0 & PC2 > 0 ~ -0.08,   # quadrant -+
      PC1 < 0 & PC2 < 0 ~ -0.08,   # quadrant --
      PC1 > 0 & PC2 < 0 ~  0.08,   # quadrant +-
      TRUE              ~  0       # on the axes
    ),
    nudge_y = case_when(
      PC2 > 0 ~  0.08,
      PC2 < 0 ~ -0.08,
      TRUE    ~  0
    )
  )

loadings_df_plot <- loadings_df |>
  mutate(
    PC1_plot = PC1 * 4,   # longer arrows, for readability only
    PC2_plot = PC2 * 4
  )

pca_circle <- ggplot() +

  # --- unit circle, solid line ---
  geom_path(
    data  = circle_df,
    aes(x = x, y = y),
    color = "grey40",
    linewidth = 0.8,
    linetype = "solid"
  ) +

  # --- quadrant axes, solid lines ---
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +

  # --- loading arrows ---
  geom_segment(
    data = loadings_df_plot,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = magnitude),
    arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.1
  ) +

  # --- trait labels ---
  geom_text_repel(
    data              = loadings_df_plot,
    aes(x = PC1, y = PC2, label = label),
    size              = 4.2,
    color             = "grey15",
    fontface          = "italic",
    box.padding       = 0.5,
    point.padding     = 0.35,
    segment.color     = "grey55",
    segment.size      = 0.35,
    max.overlaps      = Inf,
    seed              = 42,
    nudge_x           = loadings_df$nudge_x,
    nudge_y           = loadings_df$nudge_y
  ) +

  # --- color scale and theme ---
  scale_color_gradientn(
    colors   = paletteer_c("ggthemes::Temperature Diverging", 30),
    name     = "Loading\\nmagnitude",
    limits   = c(0, 0.7)
  ) +
  labs(
    x     = paste0("PC1 (", pct_pc1, "%)"),
    y     = paste0("PC2 (", pct_pc2, "%)"),
    title = "PCA correlation circle — freshwater fish morphology"
  ) +
  coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title       = element_text(size = 12),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9)
  )

print(pca_circle)

ggsave(
  filename = "figures/pca_correlation_circle_highres.png",
  plot     = pca_circle,
  width    = 10,      # inches
  height   = 10,      # inches
  units    = "in",
  dpi      = 600      # high resolution
)

# ------------------------------------------------------------------------------
# 2. Bare correlation circle (single colour, no labels)
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)

# ---- Data preparation ----
pca_trait <- readRDS("output/pca_trait.rds")
pca_obj <- pca_trait$pca_object
loadings_raw <- pca_obj$loadings[, 1:4]

# Flip PC1 so that the orientation matches the other figures
loadings_raw[, 1] <- -loadings_raw[, 1]

# Variance explained (used for the axis labels)
sdev <- pca_obj$sdev
var_explained <- sdev^2 / sum(sdev^2)
pct_pc1 <- round(var_explained[1] * 100, 1)
pct_pc2 <- round(var_explained[2] * 100, 1)

# Loadings, without labels
loadings_df <- loadings_raw |>
  as.data.frame() |>
  rownames_to_column("trait_abbr") |>
  rename(PC1 = Comp.1, PC2 = Comp.2) |>
  mutate(
    PC1_plot = PC1 * 1.6,  # arrow scaling
    PC2_plot = PC2 * 1.6
  )

# Unit circle
circle_df <- tibble(
  angle = seq(0, 2 * pi, length.out = 300),
  x = cos(angle),
  y = sin(angle)
)

# ---- Visualization ----
pca_circle_new <- ggplot() +

  # Unit circle
  geom_path(
    data = circle_df,
    aes(x = x, y = y),
    color = "grey40",
    linewidth = 0.8,
    linetype = "solid"
  ) +

  # Axes
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +

  # Loading arrows, single colour, no labels
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
    color = "black",  # same colour for every arrow
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.1
  ) +

  # Scale and theme
  labs(
    x = paste0("PC1 (", pct_pc1, "%)"),
    y = paste0("PC2 (", pct_pc2, "%)"),
    title = "PCA correlation circle — Uniform color, no labels"
  ) +
  coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title = element_text(size = 12)
  )

print(pca_circle_new)

ggsave(
  filename = "figures/pca_correlation_circle.png",
  plot = pca_circle_new,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------------------------
# 3. Heatmap of the trait-axis correlations (Figure S4)
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Trait labels ----
trait_labels <- c(
  es  = "Relative eye size",
  ep  = "Vertical eye position",
  ms  = "Relative maxillary length",
  mp  = "Oral gape position",
  elo = "Body elongation",
  wid = "Body lateral shape",
  pp  = "Pectoral fin vertical position",
  ps  = "Pectoral fin size",
  cs  = "Caudal peduncle throttling",
  svl = "Standard body length",
  bm  = "Body mass"
)

# ---- Correlations between scaled traits and PCA scores ----

# Flip PC1 so that the orientation matches the other figures
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

cor_matrix <- cor(
  pca_trait$traits_scaled,
  as.matrix(pca_trait$traits_scores),
  method = "pearson",
  use    = "pairwise.complete.obs"
)

# ---- P-values with Bonferroni correction ----
n_obs <- nrow(pca_trait$traits_scaled)
n_tests <- ncol(pca_trait$traits_scaled) * 4L  # 11 traits x 4 PC = 44

pval_matrix <- matrix(
  NA_real_,
  nrow = ncol(pca_trait$traits_scaled),
  ncol = 4,
  dimnames = list(colnames(pca_trait$traits_scaled),
                  colnames(pca_trait$traits_scores)[1:4])
)

for (tr in rownames(pval_matrix)) {
  for (pc in colnames(pval_matrix)) {
    pval_matrix[tr, pc] <- cor.test(
      pca_trait$traits_scaled[, tr],
      pca_trait$traits_scores[, pc],
      method = "pearson",
      use    = "pairwise.complete.obs"
    )$p.value
  }
}

# Bonferroni correction (adjusted p = raw p x n_tests, capped at 1)
pval_adj <- matrix(
  p.adjust(pval_matrix, method = "bonferroni"),
  nrow = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

# Significance stars, after correction
sig_vec <- dplyr::case_when(
  as.vector(pval_adj) < 0.001 ~ "***",
  as.vector(pval_adj) < 0.01  ~ "**",
  as.vector(pval_adj) < 0.05  ~ "*",
  .default = "ns"
)

sig_matrix <- matrix(
  sig_vec,
  nrow     = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

# ---- Long format ----
sig_long <- sig_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("trait") |>
  tidyr::pivot_longer(
    cols      = -trait,
    names_to  = "PC",
    values_to = "sig"
  )

cor_long <- cor_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("trait") |>
  tidyr::pivot_longer(
    cols      = -trait,
    names_to  = "PC",
    values_to = "r"
  ) |>
  dplyr::left_join(sig_long, by = c("trait", "PC")) |>
  dplyr::mutate(
    trait_label = factor(trait_labels[trait], levels = rev(trait_labels)),
    PC = factor(PC,
                levels = c("Comp.1", "Comp.2", "Comp.3", "Comp.4"),
                labels = c("PC1", "PC2", "PC3", "PC4"))
  )

# ---- Heatmap ----
p_heatmap <- ggplot(cor_long, aes(x = PC, y = trait_label, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(
      label  = sprintf("%.2f", r),
      colour = abs(r) > 0.4
    ),
    size = 3.2,
    vjust = 0.3   # slightly raised, to leave room for the star
  ) +
  geom_text(
    aes(label = sig),
    colour = "black",
    size   = 3.5,
    vjust  = -0.8  # star placed below the correlation value
  ) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey25"),
                      guide  = "none") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Pearson r"
  ) +
  scale_x_discrete(position = "top") +
  labs(title = NULL, caption = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(size = 10),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

print(p_heatmap)

ggsave(
  filename = "figures/figS4_clean.png",
  plot     = p_heatmap,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 300
)
