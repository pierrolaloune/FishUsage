# --- 13_Single_vs_MI ----------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

traitsData       <- read.table("dataPrepared/Fish/TraitFishMissing.txt", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::select(-IUCN)
traitsDataSingle <- read.table("dataPrepared/Fish/TraitFishImputed.txt", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::select(-IUCN)
selectedTraits   <- colnames(traitsDataSingle)

pcoaPhyl <- read.table("dataPrepared/Fish/pcoaPhylogenyFish.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(pcoaPhyl) <- paste0("Eigen.", 1:ncol(pcoaPhyl))
rownames(pcoaPhyl) <- gsub("Centromochlus_musaicus", "Centromochlus_musaica", rownames(pcoaPhyl))

pca_trait <- readRDS("output/pca_trait.rds")

sp_names  <- rownames(traitsData)
common_sp <- intersect(sp_names, gsub("_", " ", rownames(pcoaPhyl)))

imputation_matrix <- cbind(traitsData[, selectedTraits], pcoaPhyl[gsub(" ", "_", common_sp), ])

na_positions <- setNames(lapply(selectedTraits, function(col) which(is.na(traitsData[[col]]))), selectedTraits)

single_na_vals <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
  idx <- na_positions[[col]]
  if (length(idx) == 0) return(NULL)
  data.frame(trait = col, sp = sp_names[idx], single_val = traitsDataSingle[idx, col])
}))

ref_scores <- pca_trait$pca_object$scores[, 1:4]
mean_ref   <- attr(pca_trait$traits_scaled, "scaled:center")[selectedTraits]
sd_ref     <- attr(pca_trait$traits_scaled, "scaled:scale")[selectedTraits]
N_PC       <- 4

compute_pca_scores <- function(ximp, selectedTraits, mean_ref, sd_ref, ref_scores, sp_names, N_PC) {
  traits_scaled <- sweep(as.data.frame(ximp)[, selectedTraits, drop = FALSE], 2, mean_ref, "-")
  traits_scaled <- sweep(traits_scaled, 2, sd_ref, "/")
  scores_m <- princomp(traits_scaled)$scores[, 1:N_PC, drop = FALSE]
  rownames(scores_m) <- rownames(ximp)
  for (pc in seq_len(N_PC)) {
    if (cor(scores_m[sp_names, pc], ref_scores[, pc]) < 0) scores_m[, pc] <- -scores_m[, pc]
  }
  scores_m[sp_names, ]
}

# --- Imputations [LONG] -------------------------------------------------------

# M <- 100
# imputed_na  <- vector("list", M)
# scores_list <- vector("list", M)
#
# doParallel::registerDoParallel(cores = ncol(imputation_matrix))
#
# for (m in seq_len(M)) {
#   set.seed(100 + (m * 12))
#   ximp <- tryCatch(
#     missForest::missForest(xmis = imputation_matrix, ntree = 100, maxiter = 10,
#                            parallelize = "variables", verbose = FALSE)$ximp,
#     error = function(e) NULL
#   )
#   if (is.null(ximp)) next
#
#   imputed_na[[m]] <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
#     idx <- na_positions[[col]]
#     if (length(idx) == 0) return(NULL)
#     data.frame(imp = m, trait = col, sp = sp_names[idx], imp_value = ximp[idx, col])
#   }))
#
#   scores_list[[m]] <- tryCatch(
#     compute_pca_scores(ximp, selectedTraits, mean_ref, sd_ref, ref_scores, sp_names, N_PC),
#     error = function(e) NULL
#   )
# }
#
# doParallel::stopImplicitCluster()
#
# ok <- !sapply(imputed_na, is.null) & !sapply(scores_list, is.null)
# imputed_long <- dplyr::bind_rows(imputed_na[ok])
# scores_list  <- scores_list[ok]
#
# saveRDS(imputed_long, "output/MI_imputed_na_values.rds")
# saveRDS(scores_list,  "output/MI_scores_list.rds")

imputed_long <- readRDS("output/MI_imputed_na_values.rds")
scores_list  <- readRDS("output/MI_scores_list.rds")
M_ok         <- length(scores_list)

# --- Variability --------------------------------------------------------------

obs_ranges <- sapply(selectedTraits, function(col) diff(range(traitsData[[col]], na.rm = TRUE)))

nrmse_summary <- imputed_long %>%
  group_by(trait, sp) %>%
  summarise(sd_imp = sd(imp_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(NRMSE_sp = sd_imp / obs_ranges[trait] * 100) %>%
  group_by(trait) %>%
  summarise(
    n_missing    = n(),
    mean_NRMSE   = round(mean(NRMSE_sp, na.rm = TRUE), 3),
    median_NRMSE = round(median(NRMSE_sp, na.rm = TRUE), 3),
    sd_NRMSE     = round(sd(NRMSE_sp, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_NRMSE))

print(nrmse_summary)

# --- Single vs multiple -------------------------------------------------------

comparison_summary <- imputed_long %>%
  group_by(trait, sp) %>%
  summarise(mean_imp = mean(imp_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(single_na_vals, by = c("trait", "sp")) %>%
  mutate(diff_pct = abs(single_val - mean_imp) / obs_ranges[trait] * 100) %>%
  group_by(trait) %>%
  summarise(
    n_species       = n(),
    mean_diff_pct   = round(mean(diff_pct, na.rm = TRUE), 3),
    median_diff_pct = round(median(diff_pct, na.rm = TRUE), 3),
    sd_diff_pct     = round(sd(diff_pct, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_diff_pct))

print(comparison_summary)

# --- Procrustes [LONG] --------------------------------------------------------

# proc_df <- dplyr::bind_rows(lapply(1:(M_ok - 1), function(i) {
#   dplyr::bind_rows(lapply((i + 1):M_ok, function(j) {
#     pr <- tryCatch(
#       vegan::protest(X = scores_list[[i]], Y = scores_list[[j]], permutations = 0, symmetric = TRUE),
#       error = function(e) NULL
#     )
#     if (is.null(pr)) return(NULL)
#     data.frame(imp_i = i, imp_j = j, m2 = round(pr$ss, 6), r = round(pr$t0, 6))
#   }))
# }))
#
# write.csv(proc_df, "output/MI_procrustes_pairs.csv", row.names = FALSE)

proc_df <- read.csv("output/MI_procrustes_pairs.csv")

proc_summary <- proc_df %>%
  summarise(mean_r = mean(r, na.rm = TRUE), min_r = min(r, na.rm = TRUE), mean_m2 = mean(m2, na.rm = TRUE))

print(proc_summary)

# --- Save ---------------------------------------------------------------------

write.csv(nrmse_summary, "output/MI_NRMSE_by_trait.csv", row.names = FALSE)
write.csv(comparison_summary, "output/MI_single_vs_multiple_by_trait.csv", row.names = FALSE)

# --- Session ------------------------------------------------------------------

sessionInfo()
