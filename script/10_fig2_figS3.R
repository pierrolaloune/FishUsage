# ------------------------------------------------------------------------------
# Script: 10_fig2_figS1
# Authors: P. Bouchet, A. Toussaint
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------
# To quantify local under-representation of each usage relative to the global pool, we compute in
# each grid cell the number of species occupying that cell for ALL and for each usage. We then
# derive a bounded deficit surface Deficit(cell) = (N_ALL(cell) - N_Usage(cell)) / N_ALL(cell) in [0, 1].

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# ---- Inputs ----
tpd_trait   <- readRDS("output/TPDs_fish.rds")     # optional (kept for reproducibility)
pca_trait   <- readRDS("output/pca_trait.rds")     # must contain $uses and traits scores elsewhere
TPDs_fish   <- readRDS("output/TPD_2D.rds")        # must match object name used below
TPDc_Fish   <- readRDS("output/TPDc_Fish.rds")     # optional (kept for reproducibility)

suppressPackageStartupMessages({
  library(TPD)
  library(paletteer) # paletteer_c()
})

# ------------------------------------------------------------------------------
# Core maps and occupancy counts
# ------------------------------------------------------------------------------

# ---- Percentile-core map from TPDc ----
imageTPD_core <- function(tpd_c, thresholdPlot = 0.99) {
  
  TPDList <- tpd_c$TPDc$TPDc
  percentile_tables <- vector("list", length(TPDList))
  
  for (k in seq_along(TPDList)) {
    tmp <- cbind(index = seq_along(TPDList[[k]]), prob = TPDList[[k]])
    ord <- order(tmp[, "prob"], decreasing = TRUE)
    tmp <- tmp[ord, , drop = FALSE]
    tmp <- cbind(tmp, percentile = cumsum(tmp[, "prob"]))
    tmp <- tmp[order(tmp[, "index"]), , drop = FALSE]
    percentile_tables[[k]] <- tmp
  }
  names(percentile_tables) <- names(TPDList)
  
  xvals <- unique(tpd_c$data$evaluation_grid[, 1])
  yvals <- unique(tpd_c$data$evaluation_grid[, 2])
  
  out <- array(
    NA_real_,
    dim = c(length(xvals), length(yvals), length(TPDList)),
    dimnames = list(xvals, yvals, names(TPDList))
  )
  
  for (k in seq_along(TPDList)) {
    df <- tpd_c$data$evaluation_grid
    df$percentile <- percentile_tables[[k]][, "percentile"]
    
    for (j in seq_along(yvals)) {
      colAux <- df[df[, 2] == yvals[j], , drop = FALSE]
      out[, j, k] <- colAux$percentile
    }
    out[, , k][out[, , k] > thresholdPlot] <- NA_real_
  }
  
  out
}

# ---- Species occupancy per cell from RTPDs ----
occupancy_counts <- function(tpd_c) {
  
  RTPDs <- tpd_c$TPDc$RTPDs
  counts <- vector("list", length(RTPDs))
  
  for (k in seq_along(RTPDs)) {
    m <- RTPDs[[k]]
    m[m > 0] <- 1
    counts[[k]] <- rowSums(m)
    mode(counts[[k]]) <- "integer"
  }
  
  names(counts) <- names(RTPDs)
  counts
}

# ---- Reshape a cell vector to [Comp.1 x Comp.2] matrix ----
vec_to_mat <- function(v, eval_grid) {
  xvals <- unique(eval_grid[, 1])
  yvals <- unique(eval_grid[, 2])
  
  out <- matrix(NA_real_, nrow = length(xvals), ncol = length(yvals),
                dimnames = list(xvals, yvals))
  
  tmp <- eval_grid
  tmp$val <- v
  
  for (j in seq_along(yvals)) {
    colAux <- tmp[tmp[, 2] == yvals[j], , drop = FALSE]
    out[, j] <- colAux$val
  }
  out
}

# ------------------------------------------------------------------------------
# Compute communities once
# ------------------------------------------------------------------------------

traits_use <- pca_trait$uses
species_all <- rownames(traits_use)

usages_order <- c("All uses", "Fisheries", "Aquarium", "Game fish", "Aquaculture")
usages_order <- usages_order[usages_order %in% colnames(traits_use)]

comm <- matrix(
  0,
  nrow = 1 + length(usages_order),
  ncol = length(species_all),
  dimnames = list(c("ALL", usages_order), species_all)
)
comm["ALL", ] <- 1

for (catg in usages_order) {
  used_vec <- traits_use[species_all, catg] == 1
  comm[catg, used_vec] <- 1
}

# ---- Single TPDc call for ALL + all usages ----
TPDc_use <- TPD::TPDc(TPDs = TPDs_fish, sampUnit = comm) # long

# ---- Precompute cores and counts (once) ----
core_099  <- imageTPD_core(TPDc_use, thresholdPlot = 0.99)
core_full <- imageTPD_core(TPDc_use, thresholdPlot = 1)
counts    <- occupancy_counts(TPDc_use)

eval_grid <- TPDc_use$data$evaluation_grid
Comp.1Vec <- unique(eval_grid[, 1])
Comp.2Vec <- unique(eval_grid[, 2])

# ------------------------------------------------------------------------------
# Axis label helpers (PC1/PC2 + % if available)
# ------------------------------------------------------------------------------

get_pc_labels <- function(pca_obj) {
  # Tries multiple common slots; falls back to "PC 1" and "PC 2"
  var <- NULL
  if (!is.null(pca_obj$variance) && length(pca_obj$variance) >= 2) var <- pca_obj$variance
  if (!is.null(pca_obj$var_explained) && length(pca_obj$var_explained) >= 2) var <- pca_obj$var_explained
  if (!is.null(pca_obj$pca$variance) && length(pca_obj$pca$variance) >= 2) var <- pca_obj$pca$variance
  
  if (!is.null(var)) {
    # Accept either fractions (0-1) or percent (0-100)
    v1 <- var[1]; v2 <- var[2]
    if (max(c(v1, v2), na.rm = TRUE) <= 1) {
      v1 <- 100 * v1; v2 <- 100 * v2
    }
    xlab <- sprintf("PC 1 (%.1f%%)", v1)
    ylab <- sprintf("PC 2 (%.1f%%)", v2)
  } else {
    xlab <- "PC 1"
    ylab <- "PC 2"
  }
  
  list(xlab = xlab, ylab = ylab)
}

pc_labs <- get_pc_labels(pca_trait)

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

# ---- One panel ----
plot_one_panel_deficit <- function(catg,
                                   limX = c(-5, 5),
                                   limY = c(-7, 7),
                                   support_level = 0.99,
                                   contour_level = 0.99,
                                   ncol = 1000) {
  
  # Fixed axis labels provided by you
  xlab_pc <- "PC 1 (22.7%)"
  ylab_pc <- "PC 2 (21.1%)"
  
  imageMatALL_core <- core_099[, , "ALL"]
  
  # Build deficit from occupancy counts: (N_ALL - N_Usage) / N_ALL
  n_all <- counts[["ALL"]]
  n_use <- counts[[catg]]
  
  deficit_vec <- (n_all - n_use) / n_all
  deficit_vec[!is.finite(deficit_vec)] <- NA_real_  # handles division by 0
  deficit_vec[deficit_vec < 0] <- 0
  
  deficit_mat <- vec_to_mat(deficit_vec, eval_grid)
  
  # Mask outside global 99% core
  deficit_mat[is.na(imageMatALL_core)] <- NA_real_
  
  # Palette
  ColorRamp <- paletteer::paletteer_c("grDevices::Lajolla", ncol)
  
  # 1) Base layer (remove default Comp.1/Comp.2 labels explicitly)
  image(
    x = Comp.1Vec, y = Comp.2Vec, z = imageMatALL_core,
    xlim = limX, ylim = limY,
    col = ColorRamp,
    xaxs = "r", yaxs = "r", axes = FALSE, asp = 1,
    xlab = "", ylab = ""
  )
  
  # 2) Black polygon of global support
  cont_support <- contourLines(
    x = Comp.1Vec, y = Comp.2Vec, z = core_full[, , "ALL"],
    levels = support_level
  )
  for (fig in seq_along(cont_support)) {
    polygon(x = cont_support[[fig]]$x, y = cont_support[[fig]]$y, col = 1, border = NA)
  }
  
  # 3) Deficit overlay (again, suppress default labels)
  image(
    x = Comp.1Vec, y = Comp.2Vec, z = deficit_mat,
    xlim = limX, ylim = limY,
    col = ColorRamp,
    add = TRUE,
    xlab = "", ylab = ""
  )
  
  # 4) Global outline (solid black line)
  cont_outline <- contourLines(
    x = Comp.1Vec, y = Comp.2Vec, z = core_full[, , "ALL"],
    levels = contour_level
  )
  for (cont in cont_outline) {
    lines(cont$x, cont$y, lwd = 1.5, lty = 1, col = "black")
  }
  
  # 5) Axes + titles (PC labels only once, no overlap)
  box(which = "plot")
  axis(1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  axis(2, las = 1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  
  mtext(xlab_pc, side = 1, line = 2.2, cex = 1.0)
  mtext(ylab_pc, side = 2, line = 2.6, cex = 1.0)
  
  title(main = catg, cex.main = 1.2)
  
  invisible(ColorRamp)
}


# ---- Multi-panel figure + legend ----
plot_panel_deficit <- function(out_pdf = "figures/Figure_FunctionalShifts_clean.pdf",
                               width = 25, height = 5, pointsize = 20,
                               limX = c(-5, 5),
                               limY = c(-7, 7),
                               ncol = 1000) {
  
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  
  # Derive PNG path from PDF path
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  
  # Helper to draw all panels (avoids code duplication)
  draw_panels <- function() {
    layout(
      matrix(c(1, 2, 3, 4, 5, 6), nrow = 1, byrow = TRUE),
      widths = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.1)
    )
    
    par(mar = c(5.2, 5.2, 2.5, 0.5))
    
    for (catg in usages_order) {
      ColorRamp <- plot_one_panel_deficit(
        catg           = catg,
        limX           = limX,
        limY           = limY,
        support_level  = 0.99,
        contour_level  = 0.99,
        ncol           = ncol
      )
    }
    
    # Legend panel
    par(mar = c(5.2, 0.8, 2.5, 0.8))
    legend_image <- as.raster(matrix(rev(ColorRamp), ncol = 1))
    plot(c(0, 2), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
    
    rasterImage(legend_image, xleft = 0, ybottom = 0, xright = 1, ytop = 1)
    
    y_ticks  <- c(0, 0.25, 0.5, 0.75, 1)
    lab_ticks <- paste0(c(0, 25, 50, 75, 100), "%")
    
    segments(x0 = 1.00, x1 = 1.10, y0 = y_ticks, y1 = y_ticks, lwd = 1.2)
    text(x = 1.55, y = y_ticks, labels = lab_ticks, cex = 1)
    rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, border = "black", lwd = 1)
  }
  
  # --- PDF output ---
  pdf(out_pdf, width = width, height = height, pointsize = pointsize)
  draw_panels()
  dev.off()
  
  # --- PNG output (300 DPI, dimensions in inches → pixels) ---
  png(
    filename = out_png,
    width    = width * 300,
    height   = height * 300,
    res      = 300,
    pointsize = pointsize
  )
  draw_panels()
  dev.off()
  
  invisible(list(pdf = out_pdf, png = out_png))
}

# ---- Individual plots (square, high resolution) ----
plot_individual_deficit <- function(out_dir = "figures/individual_panels",
                                    limX = c(-5, 5),
                                    limY = c(-7, 7),
                                    ncol = 1000,
                                    png_px = 4000,
                                    png_res = 600) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (catg in usages_order) {
    # PNG (square, high resolution)
    png(
      filename = file.path(out_dir, paste0("Deficit_", gsub(" ", "_", catg), ".png")),
      width = png_px, height = png_px, res = png_res, type = "cairo"
    )
    on.exit(dev.off(), add = TRUE)
    par(mar = c(5.2, 5.2, 2.5, 0.8))
    plot_one_panel_deficit(catg = catg, limX = limX, limY = limY, ncol = ncol)
    dev.off()
    
    # PDF (square)
    pdf(
      file = file.path(out_dir, paste0("Deficit_", gsub(" ", "_", catg), ".pdf")),
      width = 7, height = 7, pointsize = 14
    )
    on.exit(dev.off(), add = TRUE)
    par(mar = c(5.2, 5.2, 2.5, 0.8))
    plot_one_panel_deficit(catg = catg, limX = limX, limY = limY, ncol = ncol)
    dev.off()
  }
  
  invisible(out_dir)
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

# 1) Multi-panel PDF
plot_panel_deficit()

# 2) Individual square high-resolution exports (PNG + PDF)
plot_individual_deficit()

# ------------------------------------------------------------------------------
# Plot Fig S1 PC3 PC4
# ------------------------------------------------------------------------------

# ---- One panel ----
plot_one_panel_deficit <- function(catg,
                                   limX = c(-5, 5),
                                   limY = c(-7, 7),
                                   support_level = 0.99,
                                   contour_level = 0.99,
                                   ncol = 1000) {
  
  # Fixed axis labels provided by you
  xlab_pc <- "PC 3 (17%)"
  ylab_pc <- "PC 4 (13.2%)"
  
  imageMatALL_core <- core_099[, , "ALL"]
  
  # Build deficit from occupancy counts: (N_ALL - N_Usage) / N_ALL
  n_all <- counts[["ALL"]]
  n_use <- counts[[catg]]
  
  deficit_vec <- 1-((n_all - n_use) / n_all)
  deficit_vec[!is.finite(deficit_vec)] <- NA_real_  # handles division by 0
  deficit_vec[deficit_vec < 0] <- 0
  deficit_vec[deficit_vec > 1] <- 1
  deficit_mat <- vec_to_mat(deficit_vec, eval_grid)
  
  # Mask outside global 99% core
  deficit_mat[is.na(imageMatALL_core)] <- NA_real_
  
  # Palette
  ColorRamp <- rev(paletteer::paletteer_c("grDevices::Lajolla", ncol))
  
  # 1) Base layer (remove default Comp.1/Comp.2 labels explicitly)
  image(
    x = Comp.3Vec, y = Comp.4Vec, z = imageMatALL_core,
    xlim = limX, ylim = limY,
    col = ColorRamp,
    xaxs = "r", yaxs = "r", axes = FALSE, asp = 1,
    xlab = "", ylab = ""
  )
  
  # 2) Black polygon of global support
  cont_support <- contourLines(
    x = Comp.3Vec, y = Comp.4Vec, z = core_full[, , "ALL"],
    levels = support_level
  )
  for (fig in seq_along(cont_support)) {
    polygon(x = cont_support[[fig]]$x, y = cont_support[[fig]]$y, col = 1, border = NA)
  }
  
  # 3) Deficit overlay (again, suppress default labels)
  image(
    x = Comp.3Vec, y = Comp.4Vec, z = deficit_mat,
    xlim = limX, ylim = limY,
    col = ColorRamp,
    add = TRUE,
    xlab = "", ylab = ""
  )
  
  # 4) Global outline (solid black line)
  cont_outline <- contourLines(
    x = Comp.3Vec, y = Comp.4Vec, z = core_full[, , "ALL"],
    levels = contour_level
  )
  for (cont in cont_outline) {
    lines(cont$x, cont$y, lwd = 1.5, lty = 1, col = "black")
  }
  
  # 5) Axes + titles (PC labels only once, no overlap)
  box(which = "plot")
  axis(1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  axis(2, las = 1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  
  mtext(xlab_pc, side = 1, line = 2.2, cex = 1.0)
  mtext(ylab_pc, side = 2, line = 2.6, cex = 1.0)
  
  title(main = catg, cex.main = 1.2)
  
  invisible(ColorRamp)
}

# ---- Multi-panel figure + legend ----
plot_panel_deficit <- function(out_pdf = "figures/Figure_FunctionalShifts_clean_PC34.pdf",
                               width = 25, height = 5, pointsize = 20,
                               limX = c(-5, 5),
                               limY = c(-7, 7),
                               ncol = 1000) {
  
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  
  pdf(out_pdf, width = width, height = height, pointsize = pointsize)
  on.exit(dev.off(), add = TRUE)
  
  layout(
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 1, byrow = TRUE),
    widths = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.1)
  )
  
  par(mar = c(5.2, 5.2, 2.5, 0.5))
  
  # Draw panels
  for (catg in usages_order) {
    ColorRamp <- plot_one_panel_deficit(
      catg = catg,
      limX = limX,
      limY = limY,
      support_level = 0.99,
      contour_level = 0.99,
      ncol = ncol
    )
  }
  
  # Legend with tick marks at 0/25/50/75/100%
  par(mar = c(5.2, 0.8, 2.5, 0.8))
  legend_image <- as.raster(matrix(ColorRamp, ncol = 1))
  plot(c(0, 2), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  
  # Color bar
  rasterImage(legend_image, xleft = 0, ybottom = 0, xright = 1, ytop = 1)
  
  # Tick positions and labels
  y_ticks <- c(1, 0.75, 0.5, 0.25, 0)
  lab_ticks <- paste0(c(100, 75, 50, 25, 0), "%")
  
  # Draw ticks (small segments) + labels
  segments(x0 = 1.00, x1 = 1.10, y0 = y_ticks, y1 = y_ticks, lwd = 1.2)
  graphics::text(x = 1.55, y = y_ticks, labels = lab_ticks, cex = 1)
  
  # Frame around legend
  rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, border = "black", lwd = 1)
  
  invisible(out_pdf)
}

# ---- Individual plots (square, high resolution) ----
plot_individual_deficit <- function(out_dir = "figures/individual_panels_34",
                                    limX = c(-5, 5),
                                    limY = c(-7, 7),
                                    ncol = 1000,
                                    png_px = 4000,
                                    png_res = 600) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (catg in usages_order) {
    # PNG (square, high resolution)
    png(
      filename = file.path(out_dir, paste0("Deficit_", gsub(" ", "_", catg), ".png")),
      width = png_px, height = png_px, res = png_res, type = "cairo"
    )
    on.exit(dev.off(), add = TRUE)
    par(mar = c(5.2, 5.2, 2.5, 0.8))
    plot_one_panel_deficit(catg = catg, limX = limX, limY = limY, ncol = ncol)
    dev.off()
    
    # PDF (square)
    pdf(
      file = file.path(out_dir, paste0("Deficit_", gsub(" ", "_", catg), ".pdf")),
      width = 7, height = 7, pointsize = 14
    )
    on.exit(dev.off(), add = TRUE)
    par(mar = c(5.2, 5.2, 2.5, 0.8))
    plot_one_panel_deficit(catg = catg, limX = limX, limY = limY, ncol = ncol)
    dev.off()
  }
  
  invisible(out_dir)
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

# 1) Multi-panel PDF
plot_panel_deficit()

# 2) Individual square high-resolution exports (PNG + PDF)
plot_individual_deficit()
