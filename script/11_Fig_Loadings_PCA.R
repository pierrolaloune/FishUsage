# ------------------------------------------------------------------------------
# Script : 11_Fig_Loadings_PCA
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads a princomp PCA object and computes trait–axis correlations by
# scaling PCA loadings with axis standard deviations. It produces both wide and
# long correlation tables, generates a publication-ready heatmap of correlations,
# and creates a PCA correlation circle (PC1–PC2) using correlation vectors, with
# clean exports to the figures directory.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# ---- Inputs ----
pca_trait <- readRDS("output/pca_trait.rds")
stopifnot(exists("pca_trait"))

stopifnot(!is.null(pca_trait$pca_object))
stopifnot(inherits(pca_trait$pca_object, "princomp"))

dims_used <- pca_trait$dimensions_used
stopifnot(is.numeric(dims_used), length(dims_used) == 1)
stopifnot(dims_used == 4)

# ------------------------------------------------------------------------------
# Processing
# ------------------------------------------------------------------------------

# ---- Flip PC1 for interpretability (keep everything consistent) ----
pca_trait$pca_object$scores[, 1]   <- -pca_trait$pca_object$scores[, 1]
pca_trait$pca_object$loadings[, 1] <- -pca_trait$pca_object$loadings[, 1]

if (!is.null(pca_trait$loadings)) {
  pca_trait$loadings[, 1] <- -pca_trait$loadings[, 1]
}

if (!is.null(pca_trait$traits_scores)) {
  if ("Comp.1" %in% names(pca_trait$traits_scores)) {
    pca_trait$traits_scores$Comp.1 <- -pca_trait$traits_scores$Comp.1
  }
  if ("PC1" %in% names(pca_trait$traits_scores)) {
    pca_trait$traits_scores$PC1 <- -pca_trait$traits_scores$PC1
  }
}

# 1) Extract loadings for the selected axes
loadings_mat <- unclass(pca_trait$loadings)[, 1:dims_used, drop = FALSE]
stopifnot(is.matrix(loadings_mat))
stopifnot(!anyNA(loadings_mat))

# 2) Extract sdev for matching axes
sdev_vec <- pca_trait$pca_object$sdev[1:dims_used]
stopifnot(length(sdev_vec) == ncol(loadings_mat))
stopifnot(!anyNA(sdev_vec))

# 3) Trait–axis correlations
cor_mat <- sweep(loadings_mat, MARGIN = 2, STATS = sdev_vec, FUN = "*")

# ---- Rename axes: Comp.X -> PCX (for outputs/figure only) ----
axis_names <- paste0("PC", seq_len(dims_used))
colnames(cor_mat) <- axis_names

# ------------------------------------------------------------------------------
# Outputs: tables
# ------------------------------------------------------------------------------

# Wide table (traits x axes)
cor_wide <- as.data.frame(cor_mat)
cor_wide$trait <- rownames(cor_wide)
cor_wide <- cor_wide[, c("trait", colnames(cor_mat)), drop = FALSE]
rownames(cor_wide) <- NULL

# Long table (trait, axis, correlation)
cor_long <- data.frame(
  trait = rep(rownames(cor_mat), times = ncol(cor_mat)),
  axis  = rep(colnames(cor_mat), each  = nrow(cor_mat)),
  corr  = as.vector(cor_mat),
  row.names = NULL
)

# Order traits by correlation on PC1
trait_order <- cor_long[cor_long$axis == "PC1", c("trait", "corr")]
trait_order <- trait_order[order(trait_order$corr, decreasing = FALSE), "trait"]
cor_long$trait <- factor(cor_long$trait, levels = trait_order)

# Ensure axis ordering is stable
cor_long$axis <- factor(cor_long$axis, levels = axis_names)

# ------------------------------------------------------------------------------
# Visualization: heatmap
# ------------------------------------------------------------------------------

p_heat <- ggplot(cor_long, aes(x = axis, y = trait, fill = corr)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 3.2, color = "black") +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "white", high = "#C53030",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Correlation"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

out_dir <- "figures/Clean/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir, "figS4_clean.pdf"),
  plot = p_heat,
  width = 7.8, height = 5.6, units = "in",
  dpi = 300
)

print(p_heat)

# ------------------------------------------------------------------------------
# Visualization: PCA correlation circle (PC1-PC2)
# ------------------------------------------------------------------------------

cor12 <- as.data.frame(cor_mat[, c("PC1", "PC2"), drop = FALSE])
cor12$trait <- rownames(cor12)
rownames(cor12) <- NULL

cor12$rho <- sqrt(cor12$PC1^2 + cor12$PC2^2)
cor12 <- cor12[order(cor12$rho, decreasing = TRUE), ]

theta <- seq(0, 2 * pi, length.out = 400)
circle_df <- data.frame(x = cos(theta), y = sin(theta))

p_circle <- ggplot() +
  geom_path(
    data = circle_df,
    aes(x = x, y = y),
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.5) +
  geom_segment(
    data = cor12,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.7,
    color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = cor12,
    aes(x = PC1, y = PC2, label = trait),
    size = 3.4,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 123
  ) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    x = "PC1 (correlation)",
    y = "PC2 (correlation)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "figS4b_pca_circle_PC1_PC2_clean.png"),
  plot = p_circle,
  width = 6.5, height = 6.5, units = "in",
  dpi = 300
)

print(p_circle)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# sessionInfo()
