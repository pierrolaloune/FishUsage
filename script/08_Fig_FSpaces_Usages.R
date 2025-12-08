# ============================================================
# Script: 08_Fig_FSpaces_Usages.R
# Purpose: Generate functional space (FS) plots for all human
#          use categories in PC1–PC2 and PC3–PC4 spaces
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

funspace_results <- readRDS("output/funspace_results.rds")
pca_trait <- readRDS("output/pca_trait.rds")
uni_clean <- readRDS("output/uni_clean.rds")

################################################################################
# SELECT SPECIES FOR ANNOTATION
################################################################################

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
  "Oreochromis andersonii"
)
pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]
pca_scores_selected <- pca_trait$traits_scores[matching_species, 1:4, drop = FALSE]
print(pca_scores_selected)

################################################################################
# PC1–PC2: PLOTS
################################################################################

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
  colors = colorRampPalette(c("white","#D2FF28","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  # xlim = x_limits,
  # ylim = y_limits
) 
 
 points(pca_scores_selected[, 1], pca_scores_selected[, 2], 
        col = "black", pch = 16, cex = 0.4)

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#5EB1BF","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  


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
  colors = colorRampPalette(c("#FFFFFF","#999999","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 


dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#63A088", "#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#D496A7","#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#8D86C9","#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  

################################################################################
# PC3–PC4: PLOTS
################################################################################

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
  colors = colorRampPalette(c("white","#D2FF28","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
) 

points(pca_scores_selected[, 1], pca_scores_selected[, 2], 
       col = "black", pch = 16, cex = 0.4)

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#5EB1BF","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  


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
  colors = colorRampPalette(c("#FFFFFF","#999999","#C84C09"))(1000),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 


dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#63A088", "#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#D496A7","#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  

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
  colors = colorRampPalette(c("#FFFFFF","#8D86C9","#C84C09"))(500),
  globalContour = TRUE, 
  globalContour.quant	= NULL,
  globalContour.lty = 3,                         
  globalContour.col = "black",
  globalContour.lwd	= 1,
  xlim = x_limits,
  ylim = y_limits
) 

dev.off()  
