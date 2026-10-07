# --- 08_Fig_FSpaces_Usages ----------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")
source("script/000_layout.R")

# --- Data ---------------------------------------------------------------------

funspace_results  <- readRDS("output/funspace_results.rds")
pca_trait         <- readRDS("output/pca_trait.rds")
phylopic_manifest <- read.csv("dataPrepared/phylopic_manifest.csv")

matching_species <- c(
  "Psephurus gladius",
  "Atractosteus spatula",
  "Anguilla anguilla",
  "Luciobarbus brachycephalus",
  "Wallago attu",
  "Salmo trutta",
  "Dermogenys pusilla",
  "Huso huso",
  "Oreochromis andersonii",
  "Hemiancistrus medians"
)

fs_uses <- tribble(
  ~use,          ~col,      ~n_col,
  "Alluses",     "#D2FF28", 1000,
  "Fisheries",   "#5EB1BF", 1000,
  "Aquaculture", "#999999", 1000,
  "Aquarium",    "#63A088",  500,
  "Gamefish",    "#D496A7",  500
)

pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]
pca_scores_selected <- pca_trait$traits_scores[matching_species, 1:4, drop = FALSE]

# --- Functional spaces --------------------------------------------------------

fs_panel <- function(use, axes, limits = TRUE) {
  u   <- fs_uses[fs_uses$use == use, ]
  lim <- lapply(paste0("Comp.", axes), function(k) range(pca_trait$traits_scores[, k]) * c(1.1, 1.1))
  panel_png(function() {
    plot(
      funspace_results[[paste0("FS_", use, "_PC", axes[1], "PC", axes[2])]],
      type = "groups", which.group = "1", quant.plot = TRUE,
      pnt = TRUE, pnt.col = rgb(0, 0, 0, 0.01),
      colors = colorRampPalette(c("#FFFFFF", u$col, "#C84C09"))(u$n_col),
      globalContour = TRUE, globalContour.quant = NULL, globalContour.lty = 3,
      globalContour.col = "black", globalContour.lwd = 1,
      xlim = if (limits) lim[[1]], ylim = if (limits) lim[[2]]
    )
    if (use == "Alluses") {
      points(pca_scores_selected[, axes[1]], pca_scores_selected[, axes[2]], col = "black", pch = 16, cex = 0.4)
    }
  }, 3000, 3000, res = 300, device = "png")
}

fs_png <- list(
  fig1  = purrr::map(setNames(fs_uses$use, fs_uses$use), ~ fs_panel(.x, c(1, 2))),
  figS1 = purrr::map(setNames(fs_uses$use, fs_uses$use), ~ fs_panel(.x, c(3, 4), limits = .x != "Alluses"))
)

# --- Correlation circles ------------------------------------------------------

circle_png <- list(
  fig1  = panel_png(plot_cor_circle(pca_trait, 1, 2), 6000, 6000, res = 600),
  figS1 = panel_png(plot_cor_circle(pca_trait, 3, 4, stretch = 1.3), 6000, 6000, res = 600)
)

# --- Silhouettes --------------------------------------------------------------

silhouettes <- fetch_phylopic_cache(phylopic_manifest, "dataPrepared/phylopic_silhouettes.rds")

# --- Save ---------------------------------------------------------------------

phylopic_manifest %>%
  mutate(url = paste0("https://www.phylopic.org/images/", uuid)) %>%
  write.csv("output/phylopic_credits.csv", row.names = FALSE)

fig_save("figures/fig1",  width = 439, height = 519, draw = draw_fs_figure("fig1",  fs_png, circle_png, silhouettes))
fig_save("figures/figS1", width = 422, height = 519, draw = draw_fs_figure("figS1", fs_png, circle_png, silhouettes))

# --- Session ------------------------------------------------------------------

sessionInfo()
