# ------------------------------------------------------------------------------
# Script : 05_Distinctiveness_IUCN
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads fish TPDs and the PCA/use object, builds a community matrix
# (uses × species), and computes functional uniqueness (Ui) and functional
# distinctiveness from Euclidean trait-space distances. It then merges uniqueness
# and distinctiveness at the species level, visualizes their relationship, runs a
# Pearson correlation test, and provides save lines for downstream reuse.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait      <- readRDS("output/TPDs_fish.rds")
pca_trait      <- readRDS("output/pca_trait.rds")
species_scores <- as.data.frame(pca_trait$traits_scores[, 1:4])
species_uses   <- as.data.frame(pca_trait$uses)

# ------------------------------------------------------------------------------
# Community matrix (uses × species)
# ------------------------------------------------------------------------------

species_uses$rownames   <- rownames(species_uses)
species_scores$rownames <- rownames(species_scores)

species_scores_uses <- merge(species_uses, species_scores, by = "rownames")
rownames(species_scores_uses) <- species_scores_uses$rownames

species_scores_uses$rownames <- NULL
species_uses$rownames        <- NULL

MatriceFish <- t(as.matrix(species_uses))

column_names <- rownames(species_uses)
MatriceFish_1 <- matrix(1, nrow = 1, ncol = length(column_names))
colnames(MatriceFish_1) <- column_names
rownames(MatriceFish_1) <- "all"

MatriceFish <- rbind(MatriceFish, MatriceFish_1)

# write.csv(MatriceFish, "dataPrepared/MatriceFish.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# Uniqueness calculation
# ------------------------------------------------------------------------------

species_scores_mat <- as.matrix(pca_trait$traits_scores[, 1:4])
dist_matrix <- as.matrix(dist(species_scores_mat, method = "euclidean"))

uni <- funrar::uniqueness(MatriceFish, dist_matrix)

# ------------------------------------------------------------------------------
# Distinctiveness calculation
# ------------------------------------------------------------------------------

dist_res <- funrar::distinctiveness(
  MatriceFish["all", , drop = FALSE],
  dist_matrix
)

dist_vec <- as.numeric(dist_res[1, ])
names(dist_vec) <- colnames(MatriceFish)

dist_df <- data.frame(
  species = colnames(dist_res),
  Dist    = dist_vec
)

# ------------------------------------------------------------------------------
# Merge uniqueness + distinctiveness
# ------------------------------------------------------------------------------

df_uni_dist <- data.frame(
  species         = uni$species,
  Ui              = uni$Ui,
  Distinctiveness = dist_vec
) %>%
  dplyr::filter(!is.na(Ui) & !is.na(Distinctiveness))

# ------------------------------------------------------------------------------
# Plot: uniqueness vs distinctiveness
# ------------------------------------------------------------------------------

plot_ui_dist <- ggplot2::ggplot(df_uni_dist, aes(x = Ui, y = Distinctiveness)) + # Figure S2
  ggplot2::geom_point(alpha = 0.5, size = 2) +
  ggplot2::geom_smooth(method = "lm", color = "orange", se = TRUE) +
  ggplot2::labs(
    x = "Ui (Functional Uniqueness)",
    y = "Functional Distinctiveness",
    title = "Relationship between Ui and Distinctiveness"
  ) +
  ggplot2::theme_minimal(base_size = 13)

# ggsave("figures/plot_ui_dist.jpg", plot_ui_dist, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# Correlation test
# ------------------------------------------------------------------------------

cor_test <- cor.test(df_uni_dist$Ui, df_uni_dist$Distinctiveness, method = "pearson")

r_value <- cor_test$estimate
p_value <- cor_test$p.value

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(uni, "dataPrepared/Fish/uni.rds")
# saveRDS(dist_df, "dataPrepared/Fish/dist.rds")
