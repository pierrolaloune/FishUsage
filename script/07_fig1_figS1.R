# ------------------------------------------------------------------------------
# Script : 07_fig1_FS_Uses
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads precomputed funspace objects, the PCA/use object, and a cleaned
# uniqueness table to select a subset of rare/unique used species for annotation.
# It then generates functional space plots for each human-use category in two PCA
# subspaces (PC1-PC2 and PC3-PC4), exporting each panel as a high-resolution PNG
# with a global contour overlay.
#
# Every panel follows the same recipe, only the funspace object and the colour
# ramp change:
#   All uses    #D2FF28      Fisheries   #5EB1BF      Aquaculture #999999
#   Aquarium    #63A088      Game fish   #D496A7      Bait        #8D86C9
# The dotted black line is the contour of the global functional space, so each
# use can be read against the same reference.
#
# Requires output/funspace_results.rds, produced by 02_FSpaces_Usages.R.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
funspace_results <- readRDS("output/funspace_results.rds")
pca_trait <- readRDS("output/pca_trait.rds")
uni_clean <- readRDS("output/uni_clean.rds")

# ------------------------------------------------------------------------------
# Select species for annotation
# ------------------------------------------------------------------------------

# The 100 most functionally unique used species, from which a handful of
# emblematic ones are picked by hand to be marked on the "All uses" panel.

uni_rare_clean <- uni_clean %>%
  filter(`All uses` == 1) %>%
  arrange(desc(Ui)) %>%
  slice_head(n = 100)

matching_species <- c(
  "Psephurus gladius",
  "Atractosteus spatula",
  "Anguilla anguilla",
  "Luciobarbus brachycephalus",
  "Wallago attu",
  "Chitala blanci",
  "Dermogenys pusilla",
  "Huso huso",
  "Oreochromis andersonii",
  "Hemiancistrus medians"
)

# ---- Flip PC1 so that large-bodied species sit on the right ----
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]
pca_scores_selected <- pca_trait$traits_scores[matching_species, 1:4, drop = FALSE]
print(pca_scores_selected)

# ------------------------------------------------------------------------------
# PC1-PC2: plots
# ------------------------------------------------------------------------------

# ---- All uses (annotated with the selected species) ----
png("figures/FS_Alluses_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Alluses_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("white", "#D2FF28", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  # xlim = x_limits,
  # ylim = y_limits
)

points(pca_scores_selected[, 1], pca_scores_selected[, 2],
       col = "black", pch = 16, cex = 0.4)

dev.off()

# ---- Fisheries ----
png("figures/FS_Fisheries_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Fisheries_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#5EB1BF", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Aquaculture ----
png("figures/FS_Aquaculture_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Aquaculture_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#999999", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Aquarium ----
png("figures/FS_Aquarium_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Aquarium_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#63A088", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Game fish ----
png("figures/FS_Gamefish_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Gamefish_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#D496A7", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Bait ----
png("figures/FS_Bait_PC1PC2.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.1"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.2"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Bait_PC1PC2,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#8D86C9", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ------------------------------------------------------------------------------
# PC3-PC4: plots
# ------------------------------------------------------------------------------

# Same six panels, projected on the third and fourth PCA axes.

# ---- All uses (annotated with the selected species) ----
png("figures/FS_Alluses_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Alluses_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("white", "#D2FF28", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
)

points(pca_scores_selected[, 1], pca_scores_selected[, 2],
       col = "black", pch = 16, cex = 0.4)

dev.off()

# ---- Fisheries ----
png("figures/FS_Fisheries_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Fisheries_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#5EB1BF", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Aquaculture ----
png("figures/FS_Aquaculture_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Aquaculture_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#999999", "#C84C09"))(1000),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Aquarium ----
png("figures/FS_Aquarium_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Aquarium_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#63A088", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Game fish ----
png("figures/FS_Gamefish_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Gamefish_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#D496A7", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ---- Bait ----
png("figures/FS_Bait_PC3PC4.png",
    width = 3000, height = 3000, res = 300)
x_limits <- range(pca_trait$traits_scores[, "Comp.3"]) * c(1.1, 1.1)
y_limits <- range(pca_trait$traits_scores[, "Comp.4"]) * c(1.1, 1.1)
plot(
  funspace_results$FS_Bait_PC3PC4,
  type = "groups",
  which.group = "1",
  quant.plot = TRUE,
  pnt = TRUE,
  pnt.col = rgb(0, 0, 0, 0.01),
  colors = colorRampPalette(c("#FFFFFF", "#8D86C9", "#C84C09"))(500),
  globalContour = TRUE,
  globalContour.quant = NULL,
  globalContour.lty = 3,
  globalContour.col = "black",
  globalContour.lwd = 1,
  xlim = x_limits,
  ylim = y_limits
)

dev.off()

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

# saveRDS(funspace_results, "output/funspace_results.rds")
