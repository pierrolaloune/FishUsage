# ------------------------------------------------------------------------------
# 12_Single_vs_MI
# P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Data import + preparation
# ------------------------------------------------------------------------------

traitsData <- read.table(
  "dataPrepared/Fish/TraitFishMissing.txt",
  header = TRUE, stringsAsFactors = FALSE
) %>% dplyr::select(-IUCN)

traitsDataSingle <- read.table(
  "dataPrepared/Fish/TraitFishImputed.txt",
  header = TRUE, stringsAsFactors = FALSE
) %>% dplyr::select(-IUCN)

selectedTraits <- colnames(traitsDataSingle)

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

imputation_matrix <- cbind(
  traitsData[, selectedTraits],
  pcoaPhyl[gsub(" ", "_", common_sp), ]
)

na_positions <- lapply(selectedTraits, function(col) which(is.na(traitsData[[col]])))
names(na_positions) <- selectedTraits

single_na_vals <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
  idx <- na_positions[[col]]
  if (length(idx) == 0) return(NULL)
  data.frame(trait = col, sp = sp_names[idx],
             single_val = traitsDataSingle[idx, col])
}))

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
# Run
# ------------------------------------------------------------------------------

M <- 100

imputed_na  <- vector("list", M)
scores_list <- vector("list", M)
n_cores <- min(ncol(imputation_matrix))
  
cat(sprintf("Starting %d imputations (parallelize = 'variables', %d cores)...\n",
            M, n_cores))
t_start <- proc.time()["elapsed"]

doParallel::registerDoParallel(cores = n_cores)

for (m in seq_len(M)) {
  
  t_m <- proc.time()["elapsed"]
  cat(sprintf(" Imputation %d / %d ...", m, M))
  
  set.seed(100 +(m * 12) )
  
  ximp <- tryCatch(
    missForest(
      xmis        = imputation_matrix,
      ntree       = 100,
      maxiter     = 10,
      parallelize = "variables",
      verbose     = FALSE
    )$ximp,
    error = function(e) {
      message(sprintf(" ERROR: %s", e$message))
      NULL
    }
  )
  
  cat(sprintf(" done in %.1f min\n", (proc.time()["elapsed"] - t_m) / 60))
  
  if (is.null(ximp)) next
  
  # Store imputed values at NA positions (trait-level NRMSE)
  imputed_na[[m]] <- dplyr::bind_rows(lapply(selectedTraits, function(col) {
    idx <- na_positions[[col]]
    if (length(idx) == 0) return(NULL)
    data.frame(imp = m, trait = col, sp = sp_names[idx],
               imp_value = ximp[idx, col])
  }))
  
  # Recompute PCA and store scores
  scores_list[[m]] <- tryCatch(
    compute_pca_scores(ximp, selectedTraits, mean_ref, sd_ref,
                       ref_scores, sp_names, N_PC),
    error = function(e) {
      message(sprintf(" PCA ERROR m = %d: %s", m, e$message))
      NULL
    }
  )
}

doParallel::stopImplicitCluster()

elapsed_total <- round(proc.time()["elapsed"] - t_start, 0)
cat(sprintf("Total: %.1f min\n", elapsed_total / 60))

ok   <- !sapply(imputed_na, is.null) & !sapply(scores_list, is.null)
M_ok <- sum(ok)
cat(sprintf("%d / %d successful\n", M_ok, M))

imputed_long <- dplyr::bind_rows(imputed_na[ok])
scores_list  <- scores_list[ok]

saveRDS(imputed_long, "output/MI_imputed_na_values.rds")
saveRDS(scores_list,  "output/MI_scores_list.rds")

imputed_long <- readRDS("output/MI_imputed_na_values.rds")
scores_list <- readRDS("output/MI_scores_list.rds")

# ------------------------------------------------------------------------------
# Inter imputation
# ------------------------------------------------------------------------------

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
# Single. vs. MI
# ------------------------------------------------------------------------------

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
# Procrustes
# ------------------------------------------------------------------------------

n_pairs <- M_ok * (M_ok - 1) / 2
cat(sprintf("\nComputing Procrustes for all %d pairs...\n", n_pairs))

proc_results <- vector("list", n_pairs)
idx <- 1

for (i in 1:(M_ok - 1)) {
  for (j in (i + 1):M_ok) {
    
    proc_ij <- tryCatch(
      vegan::protest(
        X            = scores_list[[i]],
        Y            = scores_list[[j]],
        permutations = 0,
        symmetric    = TRUE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(proc_ij)) {
      proc_results[[idx]] <- data.frame(
        imp_i = i,
        imp_j = j,
        m2    = round(proc_ij$ss, 6),
        r     = round(proc_ij$t0, 6)
      )
    }
    
    idx <- idx + 1
  }
  
  if (i %% 10 == 0) {
    cat(sprintf("  Procrustes: %d / %d first indices done\n", i, M_ok - 1))
  }
}

proc_df <- dplyr::bind_rows(proc_results)
write.csv(proc_df, "output/MI_procrustes_pairs.csv", row.names = FALSE)

cat("\n=== PROCRUSTES SUMMARY (all pairs) ===\n")
cat(sprintf("Mean r  : %.4f\n", mean(proc_df$r,  na.rm = TRUE)))
cat(sprintf("Min  r  : %.4f\n", min(proc_df$r,   na.rm = TRUE)))
cat(sprintf("Mean m2 : %.4f\n", mean(proc_df$m2, na.rm = TRUE)))

cat("\n=== ANALYSIS COMPLETE ===\n")
