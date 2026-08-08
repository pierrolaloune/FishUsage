# ------------------------------------------------------------------------------
# Script : 13_Single_vs_MI
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# Multiple imputation sensitivity analysis: 100 independent missForest runs,
# parallelized over variables.
#
# The main analyses rely on a single imputation of the missing traits. This
# script checks that this choice does not drive the results, by comparing:
#   1. Inter-imputation variability - how much an imputed value moves from one
#      run to the next (NRMSE, as a percentage of the trait range).
#   2. Single vs multiple - the gap between the single imputation used in the
#      paper and the mean of the 100 runs.
#   3. Geometric stability - Procrustes correlation between the PCA spaces
#      obtained from each pair of imputations.
#
# The two heavy steps (the 100 imputations, and the pairwise Procrustes) are
# commented out; their results are stored in output/ and reloaded immediately
# afterwards, so the script runs end to end as it is.
#
# See also 100_Imputation_SI_MI.R, which runs the same analysis and adds the
# variability of the PCA scores themselves.

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(missForest)
library(vegan)
library(doParallel)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# Data import + preparation
# ------------------------------------------------------------------------------

# ---- Traits with their original missing values ----
traitsData <- read.table(
  "dataPrepared/Fish/TraitFishMissing.txt",
  header = TRUE, stringsAsFactors = FALSE
) %>% dplyr::select(-IUCN)

# ---- Reference single imputation, the one used in the paper ----
traitsDataSingle <- read.table(
  "dataPrepared/Fish/TraitFishImputed.txt",
  header = TRUE, stringsAsFactors = FALSE
) %>% dplyr::select(-IUCN)

selectedTraits <- colnames(traitsDataSingle)

# ---- Phylogenetic coordinates, used as extra predictors ----
pcoaPhyl <- read.table(
  "dataPrepared/Fish/pcoaPhylogenyFish.txt",
  header = TRUE, stringsAsFactors = FALSE
)
colnames(pcoaPhyl) <- paste0("Eigen.", 1:ncol(pcoaPhyl))
rownames(pcoaPhyl) <- gsub("Centromochlus_musaicus",
                           "Centromochlus_musaica",
                           rownames(pcoaPhyl))

sp_names  <- rownames(traitsData)
sp_phylo  <- gsub("_", " ", rownames(pcoaPhyl))
common_sp <- intersect(sp_names, sp_phylo)

# ---- Imputation matrix: traits + phylogenetic eigenvectors ----
imputation_matrix <- cbind(
  traitsData[, selectedTraits],
  pcoaPhyl[gsub(" ", "_", common_sp), ]
)

# ---- Positions of the true missing values, per trait ----
na_positions <- lapply(selectedTraits, function(col) which(is.na(traitsData[[col]])))
names(na_positions) <- selectedTraits

# ---- Values the single imputation put at those positions ----
single_na_vals <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
  idx <- na_positions[[col]]
  if (length(idx) == 0) return(NULL)
  data.frame(trait = col, sp = sp_names[idx],
             single_val = traitsDataSingle[idx, col])
}))

# ---- Reference PCA (built on the single imputation) ----
pca_trait  <- readRDS("output/pca_trait.rds")
ref_scores <- pca_trait$pca_object$scores[, 1:4]        # species x PC1-PC4
mean_ref   <- attr(pca_trait$traits_scaled, "scaled:center")[selectedTraits]
sd_ref     <- attr(pca_trait$traits_scaled, "scaled:scale")[selectedTraits]

N_PC     <- 4
n_sp     <- length(sp_names)
n_cores  <- max(1, parallel::detectCores() - 1)

# ------------------------------------------------------------------------------
# Function
# ------------------------------------------------------------------------------

# Recomputes a PCA from one imputed table and returns the species scores.
# The sign of each axis is arbitrary in a PCA, so every axis is flipped when it
# is negatively correlated with the reference; without this the comparison
# between imputations would be meaningless.
compute_pca_scores <- function(ximp, selectedTraits, mean_ref, sd_ref,
                               ref_scores, sp_names, N_PC) {
  traits_imp    <- as.data.frame(ximp)[, selectedTraits, drop = FALSE]
  traits_scaled <- sweep(traits_imp, 2, mean_ref, "-")
  traits_scaled <- sweep(traits_scaled, 2, sd_ref,  "/")
  pca_m         <- princomp(traits_scaled)
  scores_m      <- pca_m$scores[, 1:N_PC, drop = FALSE]
  rownames(scores_m) <- rownames(ximp)
  for (pc in seq_len(N_PC)) {
    if (cor(scores_m[sp_names, pc], ref_scores[, pc]) < 0) {
      scores_m[, pc] <- -scores_m[, pc]
    }
  }
  scores_m[sp_names, ]
}

# ------------------------------------------------------------------------------
# Run: 100 independent missForest imputations  [LONG]
# ------------------------------------------------------------------------------

# LONG: each imputation is a full random-forest run (100 trees, 10 iterations)
# over the whole trait table. Count in hours, not minutes.

# M <- 100
#
# imputed_na  <- vector("list", M)
# scores_list <- vector("list", M)
# n_cores <- min(ncol(imputation_matrix))
#
# cat(sprintf("Starting %d imputations (parallelize = 'variables', %d cores)...\n",
#             M, n_cores))
# t_start <- proc.time()["elapsed"]
#
# doParallel::registerDoParallel(cores = n_cores)
#
# for (m in seq_len(M)) {
#
#   t_m <- proc.time()["elapsed"]
#   cat(sprintf(" Imputation %d / %d ...", m, M))
#
#   set.seed(100 + (m * 12))
#
#   ximp <- tryCatch(
#     missForest(
#       xmis        = imputation_matrix,
#       ntree       = 100,
#       maxiter     = 10,
#       parallelize = "variables",
#       verbose     = FALSE
#     )$ximp,
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
#   # Store imputed values at NA positions (trait-level NRMSE)
#   imputed_na[[m]] <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
#     idx <- na_positions[[col]]
#     if (length(idx) == 0) return(NULL)
#     data.frame(imp = m, trait = col, sp = sp_names[idx],
#                imp_value = ximp[idx, col])
#   }))
#
#   # Recompute PCA and store scores
#   scores_list[[m]] <- tryCatch(
#     compute_pca_scores(ximp, selectedTraits, mean_ref, sd_ref,
#                        ref_scores, sp_names, N_PC),
#     error = function(e) {
#       message(sprintf(" PCA ERROR m = %d: %s", m, e$message))
#       NULL
#     }
#   )
# }
#
# doParallel::stopImplicitCluster()
#
# elapsed_total <- round(proc.time()["elapsed"] - t_start, 0)
# cat(sprintf("Total: %.1f min\n", elapsed_total / 60))
#
# # Keep only the imputations that completed
# ok   <- !sapply(imputed_na, is.null) & !sapply(scores_list, is.null)
# M_ok <- sum(ok)
# cat(sprintf("%d / %d successful\n", M_ok, M))
#
# imputed_long <- dplyr::bind_rows(imputed_na[ok])
# scores_list  <- scores_list[ok]
#
# saveRDS(imputed_long, "output/MI_imputed_na_values.rds")
# saveRDS(scores_list,  "output/MI_scores_list.rds")

# ---- Recommended: load the saved results ----
imputed_long <- readRDS("output/MI_imputed_na_values.rds")
scores_list <- readRDS("output/MI_scores_list.rds")

# ------------------------------------------------------------------------------
# Analysis 1: inter-imputation variability
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
    n_missing    = n(),
    mean_NRMSE   = round(mean(NRMSE_sp,   na.rm = TRUE), 3),
    median_NRMSE = round(median(NRMSE_sp, na.rm = TRUE), 3),
    sd_NRMSE     = round(sd(NRMSE_sp,     na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_NRMSE))

print(nrmse_summary)
write.csv(nrmse_summary, "output/MI_NRMSE_by_trait.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Analysis 2: single vs multiple imputation
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
    n_species       = n(),
    mean_diff_pct   = round(mean(diff_pct,   na.rm = TRUE), 3),
    median_diff_pct = round(median(diff_pct, na.rm = TRUE), 3),
    sd_diff_pct     = round(sd(diff_pct,     na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff_pct))

print(comparison_summary)
write.csv(comparison_summary, "output/MI_single_vs_multiple_by_trait.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Analysis 3: Procrustes stability of the functional space  [LONG]
# ------------------------------------------------------------------------------

# A Procrustes correlation r close to 1 means the two PCA spaces are the same up
# to rotation, translation and scaling: the shape of the functional space does
# not depend on which imputation was used.

# LONG: M_ok * (M_ok - 1) / 2 pairwise comparisons, i.e. ~4,950 for 100 runs.

# n_pairs <- M_ok * (M_ok - 1) / 2
# cat(sprintf("\nComputing Procrustes for all %d pairs...\n", n_pairs))
#
# proc_results <- vector("list", n_pairs)
# idx <- 1
#
# for (i in 1:(M_ok - 1)) {
#   for (j in (i + 1):M_ok) {
#
#     proc_ij <- tryCatch(
#       vegan::protest(
#         X            = scores_list[[i]],
#         Y            = scores_list[[j]],
#         permutations = 0,
#         symmetric    = TRUE
#       ),
#       error = function(e) NULL
#     )
#
#     if (!is.null(proc_ij)) {
#       proc_results[[idx]] <- data.frame(
#         imp_i = i,
#         imp_j = j,
#         m2    = round(proc_ij$ss, 6),
#         r     = round(proc_ij$t0, 6)
#       )
#     }
#
#     idx <- idx + 1
#   }
#
#   if (i %% 10 == 0) {
#     cat(sprintf("  Procrustes: %d / %d first indices done\n", i, M_ok - 1))
#   }
# }
#
# proc_df <- dplyr::bind_rows(proc_results)
# write.csv(proc_df, "output/MI_procrustes_pairs.csv", row.names = FALSE)

# ---- Recommended: load the saved result ----
proc_df <- read_csv("output/MI_procrustes_pairs.csv")

cat("\n=== PROCRUSTES SUMMARY (all pairs) ===\n")
cat(sprintf("Mean r  : %.4f\n", mean(proc_df$r,  na.rm = TRUE)))
cat(sprintf("Min  r  : %.4f\n", min(proc_df$r,   na.rm = TRUE)))
cat(sprintf("Mean m2 : %.4f\n", mean(proc_df$m2, na.rm = TRUE)))

cat("\n=== ANALYSIS COMPLETE ===\n")
