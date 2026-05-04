# ------------------------------------------------------------------------------
# Script : 11_FigS3.R
# Author : P. Bouchet
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------
# This script produces two PCA visualisation figures:
#   1. A correlation circle (biplot) with full trait labels, quadrant lines,
#      and loading magnitude colour scale.
#   2. A simplified correlation circle with unlabelled arrows (single colour).
#   3. A heatmap of Pearson correlations between scaled traits and PC scores
#      (PC1–PC4), with Bonferroni-corrected significance annotations.
# PC1 is sign-flipped throughout so that higher scores correspond to larger
# body size (biological interpretability).
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
pca_trait <- readRDS("output/pca_trait.rds")
pca_obj   <- pca_trait$pca_object

# ------------------------------------------------------------------------------
# Shared data preparation
# ------------------------------------------------------------------------------

# ---- PC1 sign flip ----
# Invert PC1 so that higher values indicate larger body size (biological interpretability)
loadings_raw      <- pca_obj$loadings[, 1:4]
loadings_raw[, 1] <- -loadings_raw[, 1]

# ---- Variance explained ----
sdev          <- pca_obj$sdev
var_explained <- sdev^2 / sum(sdev^2)
pct_pc1       <- round(var_explained[1] * 100, 1)
pct_pc2       <- round(var_explained[2] * 100, 1)

# ---- Full descriptive trait labels ----
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

# ---- Unit circle ----
circle_df <- tibble(
  angle = seq(0, 2 * pi, length.out = 300),
  x     = cos(angle),
  y     = sin(angle)
)

# ------------------------------------------------------------------------------
# Figure 1 — Correlation circle with labels and magnitude colour scale
# ------------------------------------------------------------------------------

loadings_df <- loadings_raw |>
  as.data.frame() |>
  rownames_to_column("trait_abbr") |>
  rename(PC1 = Comp.1, PC2 = Comp.2) |>
  mutate(
    label     = dplyr::recode(trait_abbr, !!!trait_labels),
    magnitude = sqrt(PC1^2 + PC2^2),
    nudge_x = case_when(
      PC1 > 0 & PC2 > 0 ~  0.08,
      PC1 < 0 & PC2 > 0 ~ -0.08,
      PC1 < 0 & PC2 < 0 ~ -0.08,
      PC1 > 0 & PC2 < 0 ~  0.08,
      TRUE               ~  0
    ),
    nudge_y = case_when(
      PC2 > 0 ~  0.08,
      PC2 < 0 ~ -0.08,
      TRUE    ~  0
    )
  )

pca_circle <- ggplot() +
  
  geom_path(
    data      = circle_df,
    aes(x = x, y = y),
    color     = "grey40",
    linewidth = 0.8,
    linetype  = "solid"
  ) +
  
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +
  
  geom_segment(
    data      = loadings_df,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = magnitude),
    arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.1
  ) +
  
  geom_text_repel(
    data          = loadings_df,
    aes(x = PC1, y = PC2, label = label),
    size          = 4.2,
    color         = "grey15",
    fontface      = "italic",
    box.padding   = 0.5,
    point.padding = 0.35,
    segment.color = "grey55",
    segment.size  = 0.35,
    max.overlaps  = Inf,
    seed          = 42,
    nudge_x       = loadings_df$nudge_x,
    nudge_y       = loadings_df$nudge_y
  ) +
  
  scale_color_gradientn(
    colors = paletteer_c("ggthemes::Temperature Diverging", 30),
    name   = "Loading\nmagnitude",
    limits = c(0, 0.7)
  ) +
  
  labs(
    x     = paste0("PC1 (", pct_pc1, "%)"),
    y     = paste0("PC2 (", pct_pc2, "%)"),
    title = "PCA correlation circle — freshwater fish morphology"
  ) +
  coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title      = element_text(size = 12),
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9)
  )

print(pca_circle)

# ggsave(
#   filename = "figures/pca_correlation_circle_highres.png",
#   plot     = pca_circle,
#   width    = 10, height = 10, units = "in", dpi = 600
# )

# ------------------------------------------------------------------------------
# Figure 2 — Simplified correlation circle (no labels, single colour)
# ------------------------------------------------------------------------------

loadings_df_simple <- loadings_raw |>
  as.data.frame() |>
  rownames_to_column("trait_abbr") |>
  rename(PC1 = Comp.1, PC2 = Comp.2) |>
  mutate(
    PC1_plot = PC1 * 1.6,
    PC2_plot = PC2 * 1.6
  )

pca_circle_simple <- ggplot() +
  
  geom_path(
    data      = circle_df,
    aes(x = x, y = y),
    color     = "grey40",
    linewidth = 0.8
  ) +
  
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +
  
  geom_segment(
    data      = loadings_df_simple,
    aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
    color     = "black",
    arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 1.1
  ) +
  
  labs(
    x     = paste0("PC1 (", pct_pc1, "%)"),
    y     = paste0("PC2 (", pct_pc2, "%)"),
    title = " "
  ) +
  coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title = element_text(size = 12)
  )

print(pca_circle_simple)

# ggsave(
#   filename = "figures/pca_correlation_circle_simple.png",
#   plot     = pca_circle_simple,
#   width    = 10, height = 10, units = "in", dpi = 600
# )

# ------------------------------------------------------------------------------
# Figure 3 — Heatmap: Pearson correlations (traits × PC1–PC4)
# ------------------------------------------------------------------------------

# ---- PC1 sign flip on scores ----
# Invert PC1 so that higher values indicate larger body size (biological interpretability)
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

# ---- Pearson correlation matrix ----
cor_matrix <- cor(
  pca_trait$traits_scaled,
  as.matrix(pca_trait$traits_scores),
  method = "pearson",
  use    = "pairwise.complete.obs"
)

# ---- P-values with Bonferroni correction ----
n_tests <- ncol(pca_trait$traits_scaled) * 4L  # 11 traits × 4 PCs = 44

pval_matrix <- matrix(
  NA_real_,
  nrow     = ncol(pca_trait$traits_scaled),
  ncol     = 4,
  dimnames = list(
    colnames(pca_trait$traits_scaled),
    colnames(pca_trait$traits_scores)[1:4]
  )
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

pval_adj <- matrix(
  p.adjust(pval_matrix, method = "bonferroni"),
  nrow     = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

# ---- Significance encoding ----
sig_matrix <- matrix(
  dplyr::case_when(
    as.vector(pval_adj) < 0.001 ~ "***",
    as.vector(pval_adj) < 0.01  ~ "**",
    as.vector(pval_adj) < 0.05  ~ "*",
    .default                    ~ "ns"
  ),
  nrow     = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

# ---- Long format ----
sig_long <- sig_matrix |>
  as.data.frame() |>
  rownames_to_column("trait") |>
  pivot_longer(cols = -trait, names_to = "PC", values_to = "sig")

cor_long <- cor_matrix |>
  as.data.frame() |>
  rownames_to_column("trait") |>
  pivot_longer(cols = -trait, names_to = "PC", values_to = "r") |>
  left_join(sig_long, by = c("trait", "PC")) |>
  mutate(
    trait_label = factor(trait_labels[trait], levels = rev(trait_labels)),
    PC = factor(
      PC,
      levels = c("Comp.1", "Comp.2", "Comp.3", "Comp.4"),
      labels = c("PC1", "PC2", "PC3", "PC4")
    )
  )

# ---- Plot ----
p_heatmap <- ggplot(cor_long, aes(x = PC, y = trait_label, fill = r)) +
  
  geom_tile(colour = "white", linewidth = 0.5) +
  
  geom_text(
    aes(label = sprintf("%.2f", r), colour = abs(r) > 0.4),
    size  = 3.2,
    vjust = 0.3
  ) +
  
  geom_text(
    aes(label = sig),
    colour = "black",
    size   = 3.5,
    vjust  = -0.8
  ) +
  
  scale_colour_manual(
    values = c("TRUE" = "white", "FALSE" = "grey25"),
    guide  = "none"
  ) +
  
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
    axis.text.x     = element_text(face = "bold", size = 12),
    axis.text.y     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

print(p_heatmap)

# ggsave(
#   filename = "figures/heatmap_traits_PC.png",
#   plot     = p_heatmap,
#   width    = 5, height = 6, units = "in", dpi = 300
# )