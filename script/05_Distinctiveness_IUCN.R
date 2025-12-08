# ============================================================
# Script: 05_Distinctiveness_IUCN.R
# Purpose: Compute functional uniqueness (Ui) and distinctiveness,
#          assess their relationship, and link values to species
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

tpd_trait      <- readRDS("output/TPDs_fish.rds")
pca_trait      <- readRDS("output/pca_trait.rds")
species_scores <- as.data.frame(pca_trait$traits_scores[, 1:4])
species_uses   <- as.data.frame(pca_trait$uses)

################################################################################
# COMMUNITY MATRIX (USES × SPECIES)
################################################################################

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

################################################################################
# UNIQUENESS CALCULATION
################################################################################

species_scores_mat <- as.matrix(pca_trait$traits_scores[, 1:4])
dist_matrix <- as.matrix(dist(species_scores_mat, method = "euclidean"))

uni <- funrar::uniqueness(MatriceFish, dist_matrix)

################################################################################
# DISTINCTIVENESS CALCULATION
################################################################################

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

################################################################################
# MERGE UNIQUENESS + DISTINCTIVENESS
################################################################################

df_uni_dist <- data.frame(
  species         = uni$species,
  Ui              = uni$Ui,
  Distinctiveness = dist_vec
) %>%
  dplyr::filter(!is.na(Ui) & !is.na(Distinctiveness))

################################################################################
# PLOT: UNIQUENESS VS DISTINCTIVENESS
################################################################################

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

################################################################################
# CORRELATION TEST
################################################################################

cor_test <- cor.test(df_uni_dist$Ui, df_uni_dist$Distinctiveness, method = "pearson")

r_value <- cor_test$estimate
p_value <- cor_test$p.value

################################################################################
# SAVE OBJECTS
################################################################################

# saveRDS(uni, "dataPrepared/Fish/uni.rds")
# saveRDS(dist_df, "dataPrepared/Fish/dist.rds")
