# ------------------------------------------------------------------------------
# Script : test_new_funspace
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script computes, for each cell of the functional space, the percentage of
# species that are used, per use category, on the basis of a 99% TPD core.
#
# For every species, the cells that make up its 99% density core are identified.
# Each of those cells then gets +1 in the denominator (all species) and +1 in the
# numerator of every use the species belongs to. The ratio gives a map of the
# share of used species, cell by cell.
#
# Exploratory script: it is not part of the figures of the paper.
#
# Required objects, both loaded below:
#   TPDs_fish  : a TPDsp object (species-level TPDs on a 2D evaluation grid)
#   pca_traits : a list containing pca_traits$uses (binary 0/1 usage matrix)
#
# The plotting section also calls imageTPD(), defined in 000_functions.R.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# NOTE: absolute path, valid only on the original machine. Open the project from
# its own folder (or use an .Rproj) instead of relying on this line.
setwd("F:/00_GitHub/01_fishUsages/")

TPDs_fish <- readRDS("output/TPD_2D.rds")
pca_traits <- readRDS("output/pca_trait.rds")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Returns the indices of the smallest set of cells whose cumulative probability
# reaches 'core'. Sorting by decreasing density means the densest cells are taken
# first, which is what defines the core of the distribution.
#   prob : numeric vector (length = n_cells), should sum to 1 (or close).
#   core : numeric in (0,1], e.g. 0.99.
tpd_core_indices <- function(prob, core = 0.99) {
  stopifnot(is.numeric(prob), length(prob) > 0, core > 0, core <= 1)

  # Handle all-zero edge case (should not happen in a proper TPD)
  s <- sum(prob, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(integer(0))

  # Normalize defensively (in case of rounding)
  p <- prob / s

  ord <- order(p, decreasing = TRUE)
  cs  <- cumsum(p[ord])

  k <- which(cs >= core)[1]
  if (is.na(k)) k <- length(p)

  ord[seq_len(k)]
}

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

uses_df <- pca_traits$uses
grid_df <- TPDs_fish$data$evaluation_grid

# ------------------------------------------------------------------------------
# Processing
# ------------------------------------------------------------------------------

# ---- Harmonize species lists ----
species_tpd  <- names(TPDs_fish$TPDs)
species_uses <- rownames(uses_df)

species_common <- intersect(species_tpd, species_uses)

if (length(species_common) == 0) {
  stop("No overlapping species names between TPDs_fish$TPDs and pca_traits$uses.")
}

# Reorder to a common, reproducible order
species_common <- sort(species_common)

uses_df <- uses_df[species_common, , drop = FALSE]
tpd_list <- TPDs_fish$TPDs[species_common]

# Usage columns (all except none)
usage_names <- colnames(uses_df)

# ---- Prepare accumulators ----
n_cells <- nrow(grid_df)

denom_total <- integer(n_cells)
num_by_usage <- matrix(
  0L,
  nrow = n_cells,
  ncol = length(usage_names),
  dimnames = list(NULL, usage_names)
)

# ---- Main loop: core indices and counts, species by species  [LONG] ----
# LONG: one pass per species over the whole grid. The loop streams the counts
# rather than building a species-by-cell matrix, which would not fit in memory.
# Progress, elapsed time and ETA are printed roughly every 2 %.
core_level <- 0.99

n_sp <- length(species_common)

t0 <- Sys.time()
message(sprintf("Starting core %.2f loop for %d species...", core_level, n_sp))

# Print progress every X species
step <- max(1L, floor(n_sp / 50L))  # ~2% increments (about 50 messages max)

for (i in seq_along(species_common)) {
  sp <- species_common[i]

  prob <- unname(tpd_list[[sp]])
  idx  <- tpd_core_indices(prob, core = core_level)

  if (length(idx) == 0) next

  # Update denominator (all species)
  denom_total[idx] <- denom_total[idx] + 1L

  # Update numerator for each usage where the species is used (value == 1)
  used_cols <- usage_names[as.integer(uses_df[sp, ]) == 1L]

  if (length(used_cols) > 0) {
    for (u in used_cols) {
      num_by_usage[idx, u] <- num_by_usage[idx, u] + 1L
    }
  }

  # Progress message
  if (i %% step == 0L || i == 1L || i == n_sp) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- i / max(elapsed, 1e-9)
    remaining <- (n_sp - i) / max(rate, 1e-9)

    message(sprintf(
      "[%d/%d | %5.1f%%] elapsed: %s | eta: ~%s | last core cells: %d",
      i, n_sp, 100 * i / n_sp,
      format(as.POSIXct(elapsed, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
      format(as.POSIXct(remaining, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
      length(idx)
    ))
  }
}

message(sprintf("Finished. Total elapsed: %s",
                format(Sys.time() - t0)))

# ---- Build a long table ready for ggplot facets later ----
out_long <- as_tibble(grid_df) %>%
  mutate(cell_id = row_number(),
         n_total = denom_total) %>%
  bind_cols(as_tibble(num_by_usage) %>% rename_with(~ paste0("n_used__", .x))) %>%
  pivot_longer(
    cols = starts_with("n_used__"),
    names_to = "usage",
    values_to = "n_used"
  ) %>%
  mutate(
    usage = sub("^n_used__", "", usage),
    pct_used = if_else(n_total > 0, 100 * (n_used / n_total), NA_real_)
  ) %>%
  select(Comp.1, Comp.2, cell_id, usage, n_total, n_used, pct_used)

# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

# NOTE: writes to outputs/ (with an s), which is not the output/ directory used
# by every other script.
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

saveRDS(out_long, file = file.path("outputs", "tpd_core99_pct_used_by_cell.rds"))
write.csv(out_long, file = file.path("outputs", "tpd_core99_pct_used_by_cell.csv"), row.names = FALSE)

message("Done. Output saved to outputs/tpd_core99_pct_used_by_cell.{rds,csv}")

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

# Required object: out_long (from the previous step)
# Columns expected: Comp.1, Comp.2, usage, n_total, n_used, pct_used

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)
})

# Optional: limit to a subset of usages (comment out to keep all)
# usages_keep <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "All uses")
# plot_df <- out_long %>% filter(usage %in% usages_keep)

# ---- Quick look: base R image, one usage at a time ----
plot_one_usage_image <- function(df_long, usage_name, main_prefix = "% used (core99)") {
  df_u <- df_long %>% filter(usage == usage_name)

  x <- sort(unique(df_u$Comp.1))
  y <- sort(unique(df_u$Comp.2))

  # Build matrix z where rows correspond to x and columns correspond to y
  z <- matrix(NA_real_, nrow = length(x), ncol = length(y))
  ix <- match(df_u$Comp.1, x)
  iy <- match(df_u$Comp.2, y)
  z[cbind(ix, iy)] <- df_u$pct_used

  # Optional support matrix for contour (n_total > 0)
  z_sup <- matrix(0, nrow = length(x), ncol = length(y))
  z_sup[cbind(ix, iy)] <- ifelse(df_u$n_total > 0, 1, 0)

  # Color palette
  ncol <- 256
  cols <- viridisLite::viridis(ncol, option = "C")

  image(
    x = x, y = y, z = z,
    col = cols, zlim = c(0, 100),
    xlab = "Comp.1", ylab = "Comp.2",
    main = paste(main_prefix, "-", usage_name),
    asp = 1
  )

  # Support contour
  contour(
    x = x, y = y, z = z_sup,
    levels = 0.5, add = TRUE,
    drawlabels = FALSE,
    lwd = 0.6, col = "grey30"
  )

  box()
}

# Example:
plot_one_usage_image(out_long, "Aquaculture")

# ---- Publication style: same layout as plot_functional_shift_by_usage() ----
plot_pct_used_by_usage <- function(usage_name,
                                   df_long = out_long,
                                   TPDs = TPDs_fish,
                                   uses_df = pca_traits$uses,
                                   save_path = "figures/",
                                   file_prefix = "FS_pct_used",
                                   palette = "vik",
                                   reverse_palette = TRUE,
                                   zlim = c(0, 100),
                                   support_level = 0.999,
                                   add_usage_contour = FALSE,
                                   usage_contour_level = 0.99,
                                   width = 2000, height = 1600, res = 300) {

  # Packages needed for this plotting style
  stopifnot(requireNamespace("scico", quietly = TRUE))
  stopifnot(requireNamespace("glue", quietly = TRUE))
  stopifnot(requireNamespace("TPD", quietly = TRUE))

  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  message(glue::glue("Plotting usage: {usage_name}"))

  # --- Build matrix of pct_used like imageTPD does (column by column for Comp.2) ---
  df_u <- df_long[df_long$usage == usage_name, , drop = FALSE]
  if (nrow(df_u) == 0) stop("No rows in df_long for this usage_name: ", usage_name)

  comp1 <- sort(unique(df_u$Comp.1))
  comp2 <- sort(unique(df_u$Comp.2))

  mat_pct <- matrix(NA_real_, nrow = length(comp1), ncol = length(comp2))

  # Fill exactly like imageTPD: for each Comp.2, take the column and put it into mat
  for (i in seq_along(comp2)) {
    colAux <- df_u[df_u$Comp.2 == comp2[i], , drop = FALSE]
    colAux <- colAux[order(colAux$Comp.1), , drop = FALSE]
    mat_pct[, i] <- colAux$pct_used
  }

  # --- Recreate TPDc to get the same functional space contour ("ALL") ---
  species_all <- rownames(uses_df)
  if (is.null(species_all)) stop("uses_df must have rownames = species names.")

  if (!usage_name %in% colnames(uses_df)) {
    stop("usage_name not found in uses_df columns: ", usage_name)
  }

  used_vec <- uses_df[[usage_name]] == 1

  comm <- matrix(
    0,
    nrow = 2,
    ncol = length(species_all),
    dimnames = list(c("ALL", "Usage"), species_all)
  )
  comm["ALL", ] <- 1
  comm["Usage", used_vec] <- 1

  TPDc_use <- TPD::TPDc(TPDs = TPDs, sampUnit = comm)

  # imageTPD() is defined in 000_functions.R and must be loaded
  mat_all_full   <- imageTPD(TPDc_use, thresholdPlot = 1)[, , "ALL"]
  mat_use_full   <- imageTPD(TPDc_use, thresholdPlot = 1)[, , "Usage"]

  cont_funspace <- contourLines(
    x = unique(TPDc_use$data$evaluation_grid[, 1]),
    y = unique(TPDc_use$data$evaluation_grid[, 2]),
    z = mat_all_full,
    levels = support_level
  )

  cont_usage <- contourLines(
    x = unique(TPDc_use$data$evaluation_grid[, 1]),
    y = unique(TPDc_use$data$evaluation_grid[, 2]),
    z = mat_use_full,
    levels = usage_contour_level
  )

  # --- Colors: same aesthetic as the other maps (vik, reversed) ---
  ncol <- 1000
  ColorRamp <- scico::scico(n = ncol, palette = palette)
  if (reverse_palette) ColorRamp <- rev(ColorRamp)

  rampbreaks <- seq(zlim[1], zlim[2], length.out = ncol + 1)

  # --- Plot ---
  file_out <- glue::glue("{save_path}/{file_prefix}_{gsub(' ', '_', usage_name)}.jpg")

  jpeg(filename = file_out, width = width, height = height, res = res)

  image(
    x = comp1, y = comp2, z = mat_pct,
    col = ColorRamp, breaks = rampbreaks,
    xlab = "", ylab = "", axes = FALSE, asp = 1,
    zlim = zlim
  )

  # Contour of the global functional space
  for (cont in cont_funspace) {
    lines(cont$x, cont$y, lwd = 0.8, lty = 1, col = "grey30")
  }

  # Optional: contour of the usage core
  if (isTRUE(add_usage_contour)) {
    for (cont in cont_usage) {
      lines(cont$x, cont$y, lwd = 0.9, lty = 1, col = "black")
    }
  }

  box()
  dev.off()

  message(glue::glue("Saved: {file_out}"))
  invisible(file_out)
}

# Example:
plot_pct_used_by_usage("Aquaculture", add_usage_contour = FALSE)
