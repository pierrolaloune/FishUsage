# --- 05_Distinctiveness_IUCN --------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")
source("script/000_layout.R")

# --- Data ---------------------------------------------------------------------

pca_trait    <- readRDS("output/pca_trait.rds")
species_uses <- as.data.frame(pca_trait$uses)

MatriceFish <- rbind(
  t(as.matrix(species_uses)),
  all = rep(1, nrow(species_uses))
)

# --- Distinctiveness ----------------------------------------------------------

dist_matrix <- as.matrix(dist(pca_trait$traits_scores[, 1:4], method = "euclidean"))

uni <- funrar::uniqueness(MatriceFish, dist_matrix)

dist_res <- funrar::distinctiveness(MatriceFish["all", , drop = FALSE], dist_matrix)
dist_vec <- setNames(as.numeric(dist_res[1, ]), colnames(MatriceFish))

dist_df <- data.frame(species = colnames(dist_res), Dist = dist_vec)

df_uni_dist <- data.frame(
  species         = uni$species,
  Ui              = uni$Ui,
  Distinctiveness = dist_vec
) %>%
  filter(!is.na(Ui) & !is.na(Distinctiveness))

# --- Correlation --------------------------------------------------------------

cor_test <- cor.test(df_uni_dist$Ui, df_uni_dist$Distinctiveness, method = "pearson")

r_value <- cor_test$estimate
p_value <- cor_test$p.value

# --- Figure S5 ----------------------------------------------------------------

s5_scale <- figS5_box[["width"]] / 720

plot_ui_dist <- ggplot(df_uni_dist, aes(x = Ui, y = Distinctiveness)) +
  geom_point(alpha = 0.5, size = 2 * s5_scale, stroke = 0.5 * s5_scale) +
  geom_smooth(method = "lm", formula = y ~ x, color = "orange", se = TRUE) +
  labs(x = "Ui (Functional Uniqueness)", y = "Functional Distinctiveness", title = " ") +
  theme_minimal(base_size = 13 * s5_scale, base_family = fig_font)

label_r2 <- sprintf("R² = %.2f ; ", r_value)
label_p  <- if (p_value < 0.001) "< 0.001" else sprintf("= %.3f", p_value)

# --- Save ---------------------------------------------------------------------

# saveRDS(uni, "output/uni.rds")
# saveRDS(dist_df, "output/dist.rds")

fig_save("figures/figS5", width = 842, height = 596,
         draw = function() draw_figS5(plot_ui_dist, label_r2, label_p))

# --- Session ------------------------------------------------------------------

sessionInfo()
