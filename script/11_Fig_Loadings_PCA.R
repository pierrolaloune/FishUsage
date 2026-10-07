# --- 11_Fig_Loadings_PCA ------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")

pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

# --- Correlations -------------------------------------------------------------

cor_matrix <- cor(
  pca_trait$traits_scaled,
  as.matrix(pca_trait$traits_scores),
  method = "pearson",
  use    = "pairwise.complete.obs"
)

cor_long <- cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  pivot_longer(-trait, names_to = "PC", values_to = "r") %>%
  mutate(
    trait_label = factor(trait_labels[trait], levels = rev(trait_labels)),
    PC = factor(PC, levels = paste0("Comp.", 1:4), labels = paste0("PC", 1:4))
  )

# --- Figure S2 ----------------------------------------------------------------

s2_scale <- 842 / 576

p_heatmap <- ggplot(cor_long, aes(x = PC, y = trait_label, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5 * s2_scale) +
  geom_text(
    aes(label = sprintf("%.2f", r), colour = abs(r) > 0.4),
    size = 3.2 * s2_scale, vjust = 0.3, family = fig_font
  ) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey25"), guide = "none") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), name = "Pearson r"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11 * s2_scale, base_family = fig_font) +
  theme(
    axis.text.x     = element_text(face = "bold", size = 12 * s2_scale),
    axis.text.y     = element_text(size = 10 * s2_scale),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

# --- Save ---------------------------------------------------------------------

fig_save("figures/figS2", width = 842, height = 421,
         draw = function() fig_plot(p_heatmap, 0, 0, 842, 421))

# --- Session ------------------------------------------------------------------

sessionInfo()
