# ------------------------------------------------------------------------------
# Script : 100_Imputation_SI_MI
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script performs 100 independent missForest imputations to assess
# sensitivity to imputation uncertainty. It quantifies trait-level variability
# (NRMSE), deviation between single vs multiple imputation means, and geometric
# stability of the multidimensional functional space via Procrustes analysis.
#
# The 100 imputations are commented out; their results are stored in output/ and
# reloaded below, so the script runs end to end as it is. The Procrustes step has
# no saved fallback and does run.
#
# See also 13_Single_vs_MI.R, a near-identical version of this analysis.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# NOTE: absolute path, valid only on the original machine. Open the project from
# its own folder (or use an .Rproj) instead of relying on this line.
setwd("D:/PhD/00_R_code/01_fishUsages")

library(missForest)
library(vegan)
library(doParallel)
library(dplyr)

# ------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------

# ---- Trait data with missing values ----
traitsData <- read.table("dataPrepared/Fish/TraitFishMissing.txt",
                         header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-IUCN)

# ---- Reference single imputation ----
traitsDataSingle <- read.table("dataPrepared/Fish/TraitFishImputed.txt",
                               header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-IUCN)

selectedTraits <- colnames(traitsDataSingle)

# ---- Fixed phylogenetic coordinates ----
pcoaPhyl <- read.table("dataPrepared/Fish/pcoaPhylogenyFish.txt",
                       header = TRUE, stringsAsFactors = FALSE)
colnames(pcoaPhyl) <- paste0("Eigen.", 1:ncol(pcoaPhyl))
rownames(pcoaPhyl) <- gsub("Centromochlus_musaicus", "Centromochlus_musaica", rownames(pcoaPhyl))

sp_names <- rownames(traitsData)
sp_phylo <- gsub("_", " ", rownames(pcoaPhyl))
common_sp <- intersect(sp_names, sp_phylo)

# ---- Imputation matrix: traits + phylogenetic eigenvectors ----
imputation_matrix <- cbind(
  traitsData[, selectedTraits],
  pcoaPhyl[gsub(" ", "_", common_sp), ]
)

# ---- NA positions per trait (true missing values only) ----
na_positions <- lapply(selectedTraits, function(col) which(is.na(traitsData[[col]])))
names(na_positions) <- selectedTraits

# ---- Single imputation values at NA positions ----
single_na_vals <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
  idx <- na_positions[[col]]
  if (length(idx) == 0) return(NULL)
  data.frame(trait = col, sp = sp_names[idx], single_val = traitsDataSingle[idx, col])
}))

# ---- Reference PCA (single imputation, seed = 123) ----
pca_trait <- readRDS("output/pca_trait.rds")
ref_scores <- pca_trait$pca_object$scores[, 1:4]  # species x PC1-PC4
mean_ref <- attr(pca_trait$traits_scaled, "scaled:center")[selectedTraits]
sd_ref <- attr(pca_trait$traits_scaled, "scaled:scale")[selectedTraits]

N_PC <- 4
n_cores_imp <- min(ncol(imputation_matrix))

# ---- Helper function: recompute PCA scores ----
# The sign of a PCA axis is arbitrary, so each axis is flipped when it is
# negatively correlated with the reference; otherwise imputations could not be
# compared with one another.
compute_pca_scores <- function(ximp, selectedTraits, mean_ref, sd_ref, ref_scores, sp_names, N_PC) {
  traits_imp <- as.data.frame(ximp)[, selectedTraits, drop = FALSE]
  traits_scaled <- sweep(traits_imp, 2, mean_ref, "-")
  traits_scaled <- sweep(traits_scaled, 2, sd_ref, "/")
  pca_m <- princomp(traits_scaled)
  scores_m <- pca_m$scores[, 1:N_PC, drop = FALSE]
  rownames(scores_m) <- rownames(ximp)
  for (pc in seq_len(N_PC)) {
    if (cor(scores_m[sp_names, pc], ref_scores[, pc]) < 0) {
      scores_m[, pc] <- -scores_m[, pc]
    }
  }
  scores_m[sp_names, ]
}

# ------------------------------------------------------------------------------
# Process: 100 independent missForest imputations  [LONG]
# ------------------------------------------------------------------------------

# LONG: each imputation is a full random-forest run (100 trees, 10 iterations)
# over the whole trait table. Count in hours, not minutes.

# M <- 100
# imputed_na <- vector("list", M)
# scores_list <- vector("list", M)
#
# cat(sprintf("Starting %d imputations (parallelize = 'variables', %d cores)...\n", M, n_cores_imp))
# t_start <- proc.time()["elapsed"]
#
# registerDoParallel(cores = n_cores_imp)
#
# for (m in seq_len(M)) {
#   t_m <- proc.time()["elapsed"]
#   cat(sprintf(" Imputation %d/%d ...", m, M))
#
#   set.seed(100 + (m * 12))
#
#   ximp <- tryCatch(
#     missForest(xmis = imputation_matrix, ntree = 100, maxiter = 10,
#                parallelize = "variables", verbose = FALSE)$ximp,
#     error = function(e) {
#       message(sprintf(" ERROR: %s", e$message))
#       NULL
#     }
#   )
#
#   cat(sprintf(" done in %.1f min\n", (proc.time()["elapsed"] - t_m) / 60))
#
#   if (is.null(ximp)) next
#
#   # Store imputed values at NA positions
#   imputed_na[[m]] <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
#     idx <- na_positions[[col]]
#     if (length(idx) == 0) return(NULL)
#     data.frame(imp = m, trait = col, sp = sp_names[idx], imp_value = ximp[idx, col])
#   }))
#
#   # Recompute PCA scores
#   scores_list[[m]] <- tryCatch(
#     compute_pca_scores(ximp, selectedTraits, mean_ref, sd_ref, ref_scores, sp_names, N_PC),
#     error = function(e) {
#       message(sprintf(" PCA ERROR m = %d: %s", m, e$message))
#       NULL
#     }
#   )
# }
#
# stopImplicitCluster()
# elapsed_total <- round((proc.time()["elapsed"] - t_start) / 60, 1)
# cat(sprintf("Total: %.1f min\n", elapsed_total))
#
# # Filter successful imputations
# ok <- !sapply(imputed_na, is.null) & !sapply(scores_list, is.null)
# M_ok <- sum(ok)
# cat(sprintf("%d/%d successful\n", M_ok, M))
#
# imputed_long <- dplyr::bind_rows(imputed_na[ok])
# scores_list <- scores_list[ok]
#
# saveRDS(imputed_long, "output/MI_imputed_na_values.rds")
# saveRDS(scores_list, "output/MI_scores_list.rds")

# ---- Recommended: load the saved results ----
imputed_long <- readRDS("output/MI_imputed_na_values.rds")
scores_list <- readRDS("output/MI_scores_list.rds")
M_ok <- length(scores_list)

# ------------------------------------------------------------------------------
# Analysis 1: Inter-imputation variability (NRMSE) by trait
# ------------------------------------------------------------------------------

# How far apart the 100 runs place the same missing value, expressed as a
# percentage of the observed range of the trait.

obs_ranges <- sapply(selectedTraits, function(col) {
  diff(range(traitsData[[col]], na.rm = TRUE))
})

nrmse_summary <- imputed_long %>%
  group_by(trait, sp) %>%
  summarise(sd_imp = sd(imp_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(NRMSE_sp = sd_imp / obs_ranges[trait] * 100) %>%
  group_by(trait) %>%
  summarise(
    n_missing = n(),
    mean_NRMSE = round(mean(NRMSE_sp, na.rm = TRUE), 3),
    median_NRMSE = round(median(NRMSE_sp, na.rm = TRUE), 3),
    sd_NRMSE = round(sd(NRMSE_sp, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_NRMSE))

print("NRMSE inter-imputations by trait:")
print(nrmse_summary)
# write.csv(nrmse_summary, "output/MI_NRMSE_by_trait.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Analysis 2: Single vs multiple imputation deviation by trait
# ------------------------------------------------------------------------------

# Gap between the value used in the paper and the mean of the 100 runs, again as
# a percentage of the trait range.

comparison_summary <- imputed_long %>%
  group_by(trait, sp) %>%
  summarise(mean_imp = mean(imp_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(single_na_vals, by = c("trait", "sp")) %>%
  mutate(diff_pct = abs(single_val - mean_imp) / obs_ranges[trait] * 100) %>%
  group_by(trait) %>%
  summarise(
    n_species = n(),
    mean_diff_pct = round(mean(diff_pct, na.rm = TRUE), 3),
    median_diff_pct = round(median(diff_pct, na.rm = TRUE), 3),
    sd_diff_pct = round(sd(diff_pct, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff_pct))

print("Single vs multiple imputation differences by trait (%):")
print(comparison_summary)
# write.csv(comparison_summary, "output/MI_single_vs_multiple_by_trait.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Analysis 3: Multidimensional functional space stability
# ------------------------------------------------------------------------------

# ---- PCA scores array (species x axes x imputations) ----
scores_array <- array(NA_real_, dim = c(length(sp_names), N_PC, M_ok),
                      dimnames = list(sp_names, paste0("Comp.", 1:N_PC), paste0("m", seq_len(M_ok))))

for (m in seq_len(M_ok)) {
  scores_array[, , m] <- as.matrix(scores_list[[m]])
}

# ---- Inter-imputation variability on PCA scores ----
# The SD of a species position across imputations, relative to the length of the
# axis: how much a species moves in the functional space.
sd_array <- apply(scores_array, c(1, 2), sd, na.rm = TRUE)
axis_ranges_pca <- apply(ref_scores, 2, function(x) diff(range(x, na.rm = TRUE)))

pca_variability <- dplyr::bind_rows(lapply(1:N_PC, function(pc) {
  data.frame(
    Axis = paste0("Comp.", pc),
    axis_range = round(axis_ranges_pca[pc], 4),
    mean_SD = round(mean(sd_array[, pc], na.rm = TRUE), 6),
    mean_NRMSE_pct = round(mean(sd_array[, pc], na.rm = TRUE) / axis_ranges_pca[pc] * 100, 4),
    stringsAsFactors = FALSE
  )
}))

print("PCA scores inter-imputation variability:")
print(pca_variability)
# write.csv(pca_variability, "output/MI_pca_variability.csv", row.names = FALSE)

# ---- Procrustes analysis: all pairs  [LONG] ----
# A correlation r close to 1 means the two PCA spaces are the same up to
# rotation, translation and scaling.
# LONG: M_ok * (M_ok - 1) / 2 comparisons, i.e. ~4,950 for 100 imputations.
# Result saved as MI_procrustes_pairs.csv and reloaded below.

# n_pairs <- M_ok * (M_ok - 1) / 2
# cat(sprintf("\nComputing Procrustes for %d pairs...\n", n_pairs))
#
# proc_results <- vector("list", n_pairs)
# idx <- 1
# for (i in 1:(M_ok - 1)) {
#   for (j in (i + 1):M_ok) {
#     proc_ij <- tryCatch(
#       vegan::protest(X = scores_list[[i]], Y = scores_list[[j]],
#                      permutations = 0, symmetric = TRUE),
#       error = function(e) NULL
#     )
#     if (!is.null(proc_ij)) {
#       proc_results[[idx]] <- data.frame(
#         imp_i = i, imp_j = j,
#         m2 = round(proc_ij$ss, 6),
#         r = round(proc_ij$t0, 6)      )
#     }
#     idx <- idx + 1
#   }
# }
#
# proc_df <- dplyr::bind_rows(proc_results)

# ---- Recommended: load the saved result ----
proc_df <- read.csv("output/MI_procrustes_pairs.csv")

summary_stats <- proc_df %>%
  dplyr::summarise(
    n_pairs = nrow(.),
    mean_r = round(mean(r, na.rm = TRUE), 4),
    sd_r = round(sd(r, na.rm = TRUE), 4),
    min_r = round(min(r, na.rm = TRUE), 4),
    max_r = round(max(r, na.rm = TRUE), 4),
    mean_m2 = round(mean(m2, na.rm = TRUE), 6)
    )

print(summary_stats)
# write.csv(proc_df, "output/MI_procrustes_pairs.csv", row.names = FALSE)

cat(sprintf("Procrustes SUMMARY:\nMean r = %.4f (min = %.4f)\n",
            mean(proc_df$r, na.rm = TRUE), min(proc_df$r, na.rm = TRUE)))

cat("\n=== ANALYSIS COMPLETE ===\n")
