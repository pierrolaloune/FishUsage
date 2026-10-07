# --- 12_FS_Shifts_TPD ---------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")
source("script/000_layout.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")

pca_trait$traits_scores[, 1] <- -pca_trait$traits_scores[, 1]

usages_order <- c("All uses", "Fisheries", "Aquarium", "Aquaculture", "Game fish")

palette_pc12 <- function(n) grDevices::colorRampPalette(c("#0F1108", "#F5E0B7", "#A44200"), space = "Lab")(n)
palette_pc34 <- function(n) grDevices::hcl.colors(n, "Lajolla", rev = TRUE)

# --- Communities --------------------------------------------------------------

species_all <- rownames(pca_trait$uses)

comm <- matrix(0, nrow = 1 + length(usages_order), ncol = length(species_all),
               dimnames = list(c("ALL", usages_order), species_all))
comm["ALL", ] <- 1
for (catg in usages_order) comm[catg, pca_trait$uses[species_all, catg] == 1] <- 1

# --- TPD ----------------------------------------------------------------------

sd_pc12 <- sqrt(diag(ks::Hpi.diag(pca_trait$traits_scores[, c(1, 2)])))
sd_pc34 <- sqrt(diag(ks::Hpi.diag(cbind(-pca_trait$traits_scores[, 1], pca_trait$traits_scores[, 2]))))

tpd_plane <- function(axes, n_divisions, sds) {
  TPD::TPDsMean(
    species      = rownames(pca_trait$traits_scores),
    means        = pca_trait$traits_scores[, axes],
    sds          = matrix(rep(sds, nrow(pca_trait$traits_scores)), byrow = TRUE, ncol = 2),
    covar        = FALSE,
    alpha        = 0.95,
    samples      = NULL,
    trait_ranges = NULL,
    n_divisions  = n_divisions,
    tolerance    = 0.05
  )
}

maps_pc12 <- deficit_maps(tpd_plane(c(1, 2), 125, sd_pc12), comm)
maps_pc34 <- deficit_maps(tpd_plane(c(3, 4), 200, sd_pc34), comm)

# --- Deficit maps -------------------------------------------------------------

pad <- function(v) range(v) + c(-1, 1) * 0.05 * diff(range(v))

planes <- list(
  fig2 = list(maps = maps_pc12, limX = c(-5, 5), limY = c(-7, 7),
              xlab = "PC 1 (22.7%)", ylab = "PC 2 (21.1%)", palette = palette_pc12),
  figS3 = list(maps = maps_pc34, limX = pad(maps_pc34$eval_grid[, 1]), limY = pad(maps_pc34$eval_grid[, 2]),
               xlab = "PC 3 (17%)", ylab = "PC 4 (13.2%)", palette = palette_pc34)
)

deficit_png <- function(plane, catg) {
  panel_png(function() {
    par(mar = c(5.2, 5.2, 2.5, 0.8))
    draw_deficit_panel(plane$maps, catg, plane$limX, plane$limY, plane$xlab, plane$ylab, plane$palette)
  }, 4000, 4000, res = 600, device = "png", type = "cairo")
}

panel_files <- purrr::map(planes, function(pl) purrr::map(setNames(usages_order, usages_order), ~ deficit_png(pl, .x)))

row_file <- panel_png(function() {
  draw_deficit_row(maps_pc34, usages_order, planes$figS3$limX, planes$figS3$limY,
                   planes$figS3$xlab, planes$figS3$ylab, palette_pc34)
}, 5000, 1000, res = 200, device = "png", type = "cairo", pointsize = 20)

# --- Save ---------------------------------------------------------------------

fig_save("figures/fig2",  width = 398, height = 497, draw = draw_deficit_figure("fig2",  panel_files, row_file))
fig_save("figures/figS3", width = 376, height = 513, draw = draw_deficit_figure("figS3", panel_files, row_file))

# --- Session ------------------------------------------------------------------

sessionInfo()
