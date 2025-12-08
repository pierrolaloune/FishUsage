# ============================================================
# Script: 000_functions.R
# Purpose: Custom helper functions for TPD, null models,
#          FishBase scraping, SES, distinctiveness, GLMs,
#          and imputation error evaluation
# ============================================================

################################################################################
# TPD FUNCTIONS
################################################################################

TPDsMean_large <- function(species, means, sds, alpha = 0.95, samples = NULL,
                           trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the number of dimensions")
  }
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means' and 'sds'")
  }
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1) {
      type <- "One population_One species"
    } else {
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means' and 'sds'")
    }
    if (any(is.na(samples))) {
      stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1) {
      type <- "Multiple populations_One species"
    } else {
      type <- "Multiple populations_Multiple species"
    }
  }
  if (is.null(trait_ranges)) {
    trait_ranges <- rep(5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  if (is.null(n_divisions)) {
    n_divisions_choose <- c(1000, 200, 50, 25)
    n_divisions <- n_divisions_choose[dimensions]
  }
  grid_evaluate <- list()
  edge_length <- list()
  cell_volume <- 1
  for (dimens in 1:dimensions) {
    grid_evaluate[[dimens]] <- seq(
      from = trait_ranges[[dimens]][1],
      to   = trait_ranges[[dimens]][2],
      length = n_divisions
    )
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] - grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))) {
    names(evaluation_grid) <- paste0("Trait.", 1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions == 1) {
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  results <- list()
  results$data <- list(
    evaluation_grid = evaluation_grid,
    cell_volume = cell_volume,
    edge_length = edge_length,
    species = species,
    means = means,
    sds = sds,
    populations = if (is.null(samples)) NA else species_base,
    alpha = alpha,
    pop_means = list(),
    pop_sds = list(),
    pop_sigma = list(),
    dimensions = dimensions,
    type = type,
    method = "mean"
  )
  results$TPDs <- list()
  for (spi in 1:length(unique(species_base))) {
    if (spi == 1) {
      message(paste0("------- Calculating densities for ", type, " -----------\n"))
    }
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi] <- unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      multNormAux <- mvtnorm::dmvnorm(
        x = evaluation_grid,
        mean = results$data$pop_means[[spi]],
        sigma = results$data$pop_sigma[[spi]]
      )
      multNormAux <- multNormAux / sum(multNormAux)
      extract_alpha <- function(x) {
        alphaSpace_aux <- x[order(x, decreasing = TRUE)]
        greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux) > alpha)[1]]
        x[x < greater_prob] <- 0
        x <- x / sum(x)
        return(x)
      }
      if (alpha < 1) {
        multNormAux <- extract_alpha(multNormAux)
      }
      notZeroIndex <- which(multNormAux != 0)
      notZeroProb  <- multNormAux[notZeroIndex]
      results$TPDs[[spi]] <- cbind(notZeroIndex, notZeroProb)
    }
    if (dimensions == 1) stop("This function is intended for > 1 dimension")
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}

################################################################################
# PCA + TPD WRAPPER
################################################################################

computePCAandTPDs <- function(traits_data,
                              dimensions = NULL,
                              alpha = 0.95,
                              n_divisions_default = 100,
                              verbose = TRUE) {
  
  if (!requireNamespace("paran", quietly = TRUE)) stop("Package 'paran' needed.")
  if (!requireNamespace("TPD", quietly = TRUE)) stop("Package 'TPD' needed.")
  if (!requireNamespace("ks", quietly = TRUE)) stop("Package 'ks' needed.")
  if (!is.data.frame(traits_data)) stop("'traits_data' must be a data.frame.")
  
  traits_scaled <- scale(traits_data)
  
  if (is.null(dimensions)) {
    if (verbose) message("Estimating optimal number of dimensions using 'paran'...")
    paran_results <- paran::paran(traits_scaled, quietly = TRUE)
    dimensions <- sum(paran_results$Retained)
    if (dimensions == 0) dimensions <- 2
    if (verbose) message("Number of retained dimensions: ", dimensions)
  } else {
    if (verbose) message("Number of dimensions specified by user: ", dimensions)
  }
  
  pca_result <- princomp(traits_scaled)
  pca_summary <- summary(pca_result)$importance
  explained_variance <- (pca_summary["Standard deviation", 1:dimensions]^2) /
    sum(pca_summary["Standard deviation", ]^2)
  traits_scores <- as.data.frame(pca_result$scores[, 1:dimensions, drop = FALSE])
  colnames(traits_scores) <- paste0("Comp.", 1:dimensions)
  
  if (verbose) message("Computing TPDs...")
  grid_size <- ifelse(dimensions == 4, 30, n_divisions_default)
  sd_traits <- sqrt(diag(ks::Hpi.diag(traits_scores)))
  TPDs_result <- TPDsMean_large(
    species = rownames(traits_scores),
    means = traits_scores,
    sds = matrix(rep(sd_traits, nrow(traits_scores)), byrow = TRUE, ncol = dimensions),
    alpha = alpha,
    n_divisions = grid_size
  )
  
  output <- list(
    PCA = list(
      traits_scaled = traits_scaled,
      pca_object = pca_result,
      dimensions_used = dimensions,
      variance_explained = explained_variance,
      loadings = pca_result$loadings,
      traits_scores = traits_scores
    ),
    TPDs = TPDs_result
  )
  
  if (verbose) message("Analysis complete.")
  return(output)
}

################################################################################
# FISHBASE SCRAPING FUNCTIONS
################################################################################

make_fishbase_url <- function(species) {
  base_url <- "https://www.fishbase.se/summary/"
  species_url <- str_replace_all(tolower(species), " ", "-")
  glue("{base_url}{species_url}.html")
}

extract_human_uses <- function(url) {
  Sys.sleep(runif(1, 5, 8))
  page <- tryCatch(read_html(url), error = function(e) NULL)
  if (is.null(page)) {
    warning(glue("Page not found or connection failed for URL: {url}"))
    return(tibble(species_url = url, human_uses = NA_character_))
  }
  node_human_uses <- page %>%
    html_nodes(xpath = "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'human uses')]")
  if (length(node_human_uses) == 0) {
    human_uses_text <- NA_character_
  } else {
    sibling <- xml_find_first(node_human_uses[[1]], "following-sibling::*[1]")
    human_uses_text <- if (!is.na(sibling)) xml_text(sibling, trim = TRUE) else NA_character_
    if (is.na(human_uses_text) || human_uses_text == "") {
      parent <- xml_parent(node_human_uses[[1]])
      human_uses_text <- xml_text(parent, trim = TRUE)
      human_uses_text <- str_remove(human_uses_text, regex("Human uses[:]?	*", ignore_case = TRUE)) %>%
        str_trim()
    }
  }
  tibble(species_url = url, human_uses = human_uses_text)
}

classify_uses_precise <- function(text) {
  if (is.na(text) || text == "" || str_detect(text, fixed("Classification"))) {
    return(tibble(
      aquarium = "none",
      fisheries = "none",
      bait = "none",
      game_fish = "none",
      aquaculture = "none"
    ))
  }
  text_lower <- tolower(text)
  text_lower <- str_replace_all(text_lower, "[-]", " ")
  text_lower <- str_replace_all(text_lower, ":", ": ")
  text_lower <- str_replace_all(text_lower, " +", " ")
  
  classify_fisheries <- function(txt) {
    if (str_detect(txt, "fisheries: highly commercial")) return("highly")
    if (str_detect(txt, "fisheries: minor commercial|fisheries: subsistence fisheries")) return("rare")
    if (str_detect(txt, "fisheries: commercial")) return("regular")
    if (str_detect(txt, "fisheries: of potential interest|fisheries: of no interest")) return("none")
    if (str_detect(txt, "highly commercial") & !str_detect(txt, "aquarium|aquaculture|gamefish|bait")) return("highly")
    if (str_detect(txt, "minor commercial|subsistence fisheries")) return("rare")
    if (str_detect(txt, "commercial") & !str_detect(txt, "aquarium|aquaculture|gamefish|bait")) return("regular")
    if (str_detect(txt, "of potential interest|of no interest") & !str_detect(txt, "aquarium|aquaculture|gamefish|bait")) return("none")
    return("none")
  }
  
  classify_aquarium <- function(txt) {
    if (str_detect(txt, "aquarium: never/rarely")) return("none")
    if (str_detect(txt, "aquarium: show aquarium|aquarium: public aquarium")) return("rare")
    if (str_detect(txt, "aquarium: commercial")) return("highly")
    if (str_detect(txt, "aquarium: potential")) return("none")
    if (str_detect(txt, "never/rarely") & !str_detect(txt, "fisheries|aquaculture|gamefish|bait")) return("none")
    if (str_detect(txt, "show aquarium|public aquarium")) return("rare")
    if (str_detect(txt, "commercial") & !str_detect(txt, "fisheries|aquaculture|gamefish|bait")) return("highly")
    if (str_detect(txt, "potential") & !str_detect(txt, "fisheries|aquaculture|gamefish|bait")) return("none")
    return("none")
  }
  
  classify_aquaculture <- function(txt) {
    if (str_detect(txt, "aquaculture: commercial")) return("highly")
    if (str_detect(txt, "aquaculture: experimental")) return("rare")
    if (str_detect(txt, "aquaculture: likely future use|aquaculture: never/rarely")) return("none")
    if (str_detect(txt, "commercial") & !str_detect(txt, "fisheries|aquarium|gamefish|bait")) return("highly")
    if (str_detect(txt, "experimental")) return("rare")
    if (str_detect(txt, "likely future use|never/rarely") & !str_detect(txt, "fisheries|aquarium|gamefish|bait")) return("none")
    return("none")
  }
  
  classify_game_fish <- function(txt) {
    if (str_detect(txt, "gamefish: yes")) return("highly")
    if (str_detect(txt, "gamefish: no")) return("none")
    if (str_detect(txt, "\byes\b")) return("highly")
    if (str_detect(txt, "\bno\b")) return("none")
    return("none")
  }
  
  classify_bait <- function(txt) {
    if (str_detect(txt, "bait: usually")) return("highly")
    if (str_detect(txt, "bait: occasionally")) return("regular")
    if (str_detect(txt, "bait: never/rarely")) return("none")
    if (str_detect(txt, "usually") & !str_detect(txt, "fisheries|aquarium|gamefish|aquaculture")) return("highly")
    if (str_detect(txt, "occasionally") & !str_detect(txt, "fisheries|aquarium|gamefish|aquaculture")) return("regular")
    if (str_detect(txt, "never/rarely") & !str_detect(txt, "fisheries|aquarium|gamefish|aquaculture")) return("none")
    return("none")
  }
  
  tibble(
    aquarium   = classify_aquarium(text_lower),
    fisheries  = classify_fisheries(text_lower),
    bait       = classify_bait(text_lower),
    game_fish  = classify_game_fish(text_lower),
    aquaculture = classify_aquaculture(text_lower)
  )
}

################################################################################
# FUNCTIONAL DIVERSITY METRICS (TPDc, FRic, Dissimilarity)
################################################################################

TPDc_large <- function(TPDs, sampUnit) {
  sampUnit <- as.matrix(sampUnit)
  if (is.null(colnames(sampUnit)) | any(is.na(colnames(sampUnit)))) {
    stop("colnames(sampUnit) must contain the names of the species; NA values are not allowed")
  }
  if (is.null(rownames(sampUnit)) | any(is.na(rownames(sampUnit)))) {
    stop("rownames(sampUnit) must contain the names of the sampling units; NA values are not allowed")
  }
  if (class(TPDs) != "TPDsp") {
    stop("TPDs must be an object of class 'TPDsp', created with the function 'TPDs'")
  }
  species <- samples <- abundances <- numeric()
  for (i in 1:nrow(sampUnit)) {
    samples    <- c(samples, rep(rownames(sampUnit)[i], ncol(sampUnit)))
    species    <- c(species, colnames(sampUnit))
    abundances <- c(abundances, sampUnit[i, ])
  }
  nonZero    <- which(abundances > 0)
  samples    <- samples[nonZero]
  species    <- species[nonZero]
  abundances <- abundances[nonZero]
  results <- list()
  results$data <- TPDs$data
  results$data$sampUnit <- sampUnit
  type <- results$data$type
  if (type == "Multiple populations_One species" |
      type == "Multiple populations_Multiple species") {
    species_base <- paste(species, samples, sep = ".")
    if (!all(unique(species_base) %in% unique(results$data$populations))) {
      non_found_pops <- which(unique(species_base) %in% unique(results$data$populations) == 0)
      stop(
        "All the population TPDs must be present in 'TPDs'. Not present:\n",
        paste(species_base[non_found_pops], collapse = " / ")
      )
    }
  }
  if (type == "One population_One species" |
      type == "One population_Multiple species") {
    species_base <- species
    if (!all(unique(species_base) %in% unique(results$data$species))) {
      non_found_sps <- which(unique(species_base) %in% unique(results$data$species) == 0)
      stop(
        "All the species TPDs must be present in 'TPDs'. Not present:\n",
        paste(species_base[non_found_sps], collapse = " / ")
      )
    }
  }
  results$TPDc <- list()
  results$TPDc$species <- list()
  results$TPDc$abundances <- list()
  results$TPDc$speciesPerCell <- list()
  results$TPDc$TPDc <- list()
  
  for (samp in 1:length(unique(samples))) {
    selected_rows  <- which(samples == unique(samples)[samp])
    species_aux    <- species_base[selected_rows]
    abundances_aux <- abundances[selected_rows] / sum(abundances[selected_rows])
    RTPDsAux <- rep(0, nrow(results$data$evaluation_grid))
    TPDs_aux <- TPDs$TPDs[names(TPDs$TPDs) %in% species_aux]
    cellsOcc <- numeric()
    for (sp in 1:length(TPDs_aux)) {
      selected_name <- which(names(TPDs_aux) == species_aux[sp])
      cellsToFill   <- TPDs_aux[[selected_name]][, "notZeroIndex"]
      cellsOcc      <- c(cellsOcc, cellsToFill)
      probsToFill   <- TPDs_aux[[selected_name]][, "notZeroProb"] * abundances_aux[sp]
      RTPDsAux[cellsToFill] <- RTPDsAux[cellsToFill] + probsToFill
    }
    TPDc_aux   <- RTPDsAux
    notZeroIndex <- which(TPDc_aux != 0)
    notZeroProb  <- TPDc_aux[notZeroIndex]
    results$TPDc$TPDc[[samp]]          <- cbind(notZeroIndex, notZeroProb)
    results$TPDc$species[[samp]]       <- species_aux
    results$TPDc$abundances[[samp]]    <- abundances_aux
    results$TPDc$speciesPerCell[[samp]] <- table(cellsOcc)
    names(results$TPDc$TPDc)[samp] <-
      names(results$TPDc$species)[samp] <-
      names(results$TPDc$abundances)[samp] <-
      names(results$TPDc$speciesPerCell)[samp] <- unique(samples)[samp]
  }
  class(results) <- "TPDcomm"
  return(results)
}

Calc_FRich <- function(TPDc_Fish) {
  results_FR <- numeric()
  if (class(TPDc_Fish) == "TPDcomm") {
    TPD        <- TPDc_Fish$TPDc$TPDc
    names_aux  <- names(TPDc_Fish$TPDc$TPDc)
    cell_volume <- TPDc_Fish$data$cell_volume
  }
  if (class(TPDc_Fish) == "TPDsp") {
    TPD        <- TPDc_Fish$TPDs
    names_aux  <- names(TPDc_Fish$TPDs)
    cell_volume <- TPDc_Fish$data$cell_volume
  }
  for (i in 1:length(TPD)) {
    TPD_aux <- TPD[[i]]
    TPD_aux[TPD_aux > 0] <- cell_volume
    results_FR[i] <- sum(TPD_aux)
  }
  names(results_FR) <- names_aux
  return(results_FR)
}

dissim_large <- function(x = NULL) {
  if (class(x) == "TPDcomm") {
    TPDType <- "Communities"
    TPDc    <- x
  } else if (class(x) == "TPDsp") {
    TPDType <- "Populations"
    TPDs    <- x
  } else {
    stop("x must be an object of class TPDcomm or TPDsp")
  }
  results <- list()
  Calc_dissim <- function(x) {
    results_samp <- list()
    if (TPDType == "Communities") {
      TPD       <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
    }
    if (TPDType == "Populations") {
      TPD       <- x$TPDs
      names_aux <- names(x$TPDs)
    }
    results_samp$dissimilarity <- matrix(
      NA, ncol = length(TPD), nrow = length(TPD),
      dimnames = list(names_aux, names_aux)
    )
    results_samp$P_shared <- matrix(
      NA, ncol = length(TPD), nrow = length(TPD),
      dimnames = list(names_aux, names_aux)
    )
    results_samp$P_non_shared <- matrix(
      NA, ncol = length(TPD), nrow = length(TPD),
      dimnames = list(names_aux, names_aux)
    )
    for (i in 1:length(TPD)) {
      TPD_i <- TPD[[i]]
      for (j in 1:length(TPD)) {
        if (i > j) {
          TPD_j <- TPD[[j]]
          commonTPD <- rbind(TPD_i, TPD_j)
          duplicatedCells <- names(which(table(commonTPD[, "notZeroIndex"]) == 2))
          doubleTPD <- commonTPD[which(commonTPD[, "notZeroIndex"] %in% duplicatedCells), ]
          O_aux <- sum(tapply(doubleTPD[, "notZeroProb"], doubleTPD[, "notZeroIndex"], min))
          A_aux <- sum(tapply(doubleTPD[, "notZeroProb"], doubleTPD[, "notZeroIndex"], max)) - O_aux
          only_in_i_aux <- which(TPD_i[, "notZeroIndex"] %in%
                                   setdiff(TPD_i[, "notZeroIndex"], TPD_j[, "notZeroIndex"]))
          B_aux <- sum(TPD_i[only_in_i_aux, "notZeroProb"])
          only_in_j_aux <- which(TPD_j[, "notZeroIndex"] %in%
                                   setdiff(TPD_j[, "notZeroIndex"], TPD_i[, "notZeroIndex"]))
          C_aux <- sum(TPD_j[only_in_j_aux, "notZeroProb"])
          results_samp$dissimilarity[i, j] <- results_samp$dissimilarity[j, i] <- 1 - O_aux
          if (results_samp$dissimilarity[j, i] == 0) {
            results_samp$P_non_shared[i, j] <- NA
            results_samp$P_non_shared[j, i] <- NA
            results_samp$P_shared[i, j]     <- NA
            results_samp$P_shared[j, i]     <- NA
          } else {
            results_samp$P_non_shared[i, j] <- results_samp$P_non_shared[j, i] <-
              (2 * min(B_aux, C_aux)) / (A_aux + 2 * min(B_aux, C_aux))
            results_samp$P_shared[i, j] <- results_samp$P_shared[j, i] <-
              1 - results_samp$P_non_shared[i, j]
          }
        }
        if (i == j) {
          results_samp$dissimilarity[i, j] <- 0
        }
      }
    }
    return(results_samp)
  }
  if (TPDType == "Communities") {
    message("Computing dissimilarities between ", length(TPDc$TPDc$TPDc), " communities. This may take a while.")
    results$communities <- Calc_dissim(TPDc)
  }
  if (TPDType == "Populations") {
    message("Computing dissimilarities between ", length(TPDs$TPDs), " populations. This may take a while.")
    results$populations <- Calc_dissim(TPDs)
  }
  class(results) <- "OverlapDiss"
  return(results)
}

################################################################################
# NULL MODELS (FRic)
################################################################################

randomize_matrix <- function(original_matrix) {
  randomized <- t(apply(original_matrix, 1, function(row) sample(row)))
  colnames(randomized) <- colnames(original_matrix)
  rownames(randomized) <- rownames(original_matrix)
  return(randomized)
}

simulate_FRic_null <- function(n_iter, original_matrix, TPDs_object) {
  fric_simulations <- matrix(NA, nrow = n_iter, ncol = nrow(original_matrix))
  colnames(fric_simulations) <- rownames(original_matrix)
  for (i in 1:n_iter) {
    if (i %% 10 == 0) message("Running simulation ", i, " / ", n_iter)
    randomized_matrix <- randomize_matrix(original_matrix)
    TPDc_rand <- TPDc_large(TPDs = TPDs_object, sampUnit = randomized_matrix)
    fric_rand <- Calc_FRich(TPDc_rand)
    fric_simulations[i, ] <- fric_rand
  }
  message("Simulation complete.")
  fric_df <- as.data.frame(fric_simulations)
  fric_df$iteration <- 1:n_iter
  fric_df_long <- tidyr::pivot_longer(
    fric_df,
    cols = -iteration,
    names_to = "Usage",
    values_to = "FRic_sim"
  )
  return(fric_df_long)
}

calc_FRic_by_threat <- function(MatriceFish, TPDsp, threatsp, nrep = 999) {
  usages <- rownames(MatriceFish)
  threat_categories <- names(threatsp)
  results_list <- list()
  
  message("Starting FRic calculation by usage and threat category...")
  
  for (usage in usages) {
    message(paste0("Processing usage: ", usage))
    species_in_use <- colnames(MatriceFish)[which(MatriceFish[usage, ] == 1)]
    
    for (cat in threat_categories) {
      message(paste0("  Threat category: ", cat))
      cat_species <- threatsp[[cat]]
      
      to_remove_obs <- intersect(species_in_use, cat_species)
      n_remove <- length(to_remove_obs)
      
      if (n_remove == 0) {
        message(paste("  No species to remove for", usage, "/", cat))
        next
      }
      
      message(paste0("  Removing ", n_remove, " species for observed FRic..."))
      mat_obs <- MatriceFish[usage, , drop = FALSE]
      mat_obs[, to_remove_obs] <- 0
      TPDc_obs <- TPDc_large(TPDsp, sampUnit = mat_obs)
      FRic_obs <- Calc_FRich(TPDc_obs)[1]
      message(paste0("  Observed FRic calculated: ", round(FRic_obs, 4)))
      
      null_FRic <- numeric(nrep)
      species_in_threat <- intersect(colnames(MatriceFish), cat_species)
      message(paste0("  Launching ", nrep, " random draws in category ", cat, " (", length(species_in_threat), " possible species)..."))
      
      for (r in 1:nrep) {
        set.seed(r + 1000)
        sampled_sp <- sample(species_in_threat, n_remove)
        mat_null <- MatriceFish[usage, , drop = FALSE]
        mat_null[, sampled_sp] <- 0
        TPDc_null <- TPDc_large(TPDsp, sampUnit = mat_null)
        null_FRic[r] <- Calc_FRich(TPDc_null)[1]
        message(paste0("    Simulation ", r, " → FRic = ", round(null_FRic[r], 4)))
      }
      
      res <- data.frame(
        usage = usage,
        threat_category = cat,
        FRic_obs = FRic_obs
      )
      res[paste0("FRic_null_", 1:nrep)] <- null_FRic
      results_list[[paste(usage, cat, sep = "_")]] <- res
      message("  Results saved for this combination.\n")
    }
  }
  
  message("FRic calculation completed for all combinations.\n")
  return(do.call(rbind, results_list))
}

################################################################################
# SES (STANDARDIZED EFFECT SIZE) HELPERS
################################################################################

sesandpvalue <- function(obs, rand, nreps, probs = c(0.025, 0.975), rnd = 2) {
  if (length(rand) < 2 || all(rand == rand[1])) {
    SES <- NA
  } else {
    SES <- (obs - mean(rand)) / sd(rand)
  }
  pValsSES <- rank(c(obs, rand), ties.method = "random")[1] / (length(rand) + 1)
  results <- round(
    c(obs, SES, mean(rand), quantile(rand, prob = probs, na.rm = TRUE), pValsSES, nreps),
    rnd
  )
  names(results) <- c("Observed", "SES", "MeanRd", "CI025Rd", "CI975Rd", "Pval", "Nreps")
  return(results)
}

get_SES <- function(obs_df, sim_df, probs = c(0.025, 0.975), rnd = 2) {
  results_list <- lapply(seq_len(nrow(obs_df)), function(i) {
    usage_i <- obs_df$Use[i]
    obs_i   <- obs_df$FRich[i]
    rand_i  <- sim_df$FRic_sim[sim_df$Usage == usage_i]
    sesandpvalue(obs = obs_i, rand = rand_i, nreps = length(rand_i), probs = probs, rnd = rnd)
  })
  results_df <- as.data.frame(do.call(rbind, results_list))
  results_df$Usage <- obs_df$Use
  results_df <- dplyr::relocate(results_df, Usage)
  return(results_df)
}

plot_SES_histograms <- function(sim_df, obs_df) {
  library(ggplot2)
  library(dplyr)
  
  obs_df <- obs_df %>% rename(Usage = Use)
  sim_df <- sim_df %>% filter(Usage %in% obs_df$Usage)
  
  p <- ggplot(sim_df, aes(x = FRic_sim)) +
    geom_histogram(bins = 50, fill = "#69b3a2", alpha = 0.6, color = "grey40") +
    geom_vline(data = obs_df, aes(xintercept = FRich), color = "red", linewidth = 1) +
    facet_wrap(~Usage, scales = "free") +
    labs(
      x = "Simulated FRic", y = "Frequency",
      title = "Distribution of simulated FRic per usage",
      subtitle = "Red line = observed value"
    ) +
    theme_minimal()
  
  print(p)
}

generate_null_means <- function(pca_trait, MatriceFish, nb_simulations = 999) {
  pca_axes <- grep("^Comp\\.", colnames(pca_trait$traits_scores), value = TRUE)
  common_species <- intersect(rownames(pca_trait$traits_scores), colnames(MatriceFish))
  pca_scores <- pca_trait$traits_scores[common_species, pca_axes, drop = FALSE]
  MatriceFish <- MatriceFish[, common_species, drop = FALSE]
  result_list <- list()
  
  for (usage in rownames(MatriceFish)) {
    cat("Processing usage:", usage, "\n")
    usage_vec <- unlist(MatriceFish[usage, ])
    species_in_use <- names(usage_vec[usage_vec == 1])
    nb_species <- length(species_in_use)
    
    if (nb_species == 0) {
      warning(paste("No species associated with usage:", usage))
      next
    }
    
    observed_mean <- colMeans(pca_scores[species_in_use, , drop = FALSE])
    simulated_means <- matrix(NA, nrow = nb_simulations, ncol = length(pca_axes))
    colnames(simulated_means) <- pca_axes
    
    for (i in seq_len(nb_simulations)) {
      if (i %% 100 == 0) cat("  Simulation", i, "/", nb_simulations, "\n")
      randomized_matrix <- randomize_matrix(MatriceFish)
      usage_random_vec <- unlist(randomized_matrix[usage, ])
      species_sampled <- names(usage_random_vec[usage_random_vec == 1])
      
      if (length(species_sampled) > 0) {
        simulated_means[i, ] <- colMeans(pca_scores[species_sampled, , drop = FALSE])
      } else {
        simulated_means[i, ] <- NA
      }
    }
    
    result_list[[usage]] <- list(
      observed = observed_mean,
      simulated = simulated_means
    )
  }
  
  return(result_list)
}

get_SES_from_PCA_results <- function(results_list, probs = c(0.025, 0.975), rnd = 10) {
  output <- list()
  for (usage in names(results_list)) {
    obs_vec <- results_list[[usage]]$observed
    sim_mat <- results_list[[usage]]$simulated
    for (comp in names(obs_vec)) {
      obs_val   <- obs_vec[comp]
      rand_vals <- sim_mat[, comp]
      res <- sesandpvalue(
        obs = obs_val,
        rand = rand_vals,
        nreps = length(rand_vals),
        probs = probs,
        rnd = rnd
      )
      output[[paste(usage, comp, sep = "_")]] <- c(Usage = usage, Component = comp, res)
    }
  }
  df_out <- do.call(rbind, output)
  df_out <- as.data.frame(df_out, stringsAsFactors = FALSE)
  num_cols <- setdiff(colnames(df_out), c("Usage", "Component"))
  df_out[num_cols] <- lapply(df_out[num_cols], as.numeric)
  return(df_out)
}

calc_SES_table <- function(df, obs_col = "FRic_obs", null_prefix = "FRic_null_") {
  null_cols <- grep(paste0("^", null_prefix), names(df), value = TRUE)
  sesandpvalue_local <- function(obs, rand, nreps, probs = c(0.025, 0.975), rnd = 4) {
    SES <- (obs - mean(rand)) / sd(rand)
    pValsSES <- rank(c(obs, rand))[1] / (length(rand) + 1)
    results <- round(
      c(obs, SES, mean(rand), quantile(rand, prob = probs), pValsSES, nreps),
      rnd
    )
    names(results) <- c("Observed", "SES", "MeanRd", "CI025Rd", "CI975Rd", "Pval", "Nreps")
    return(results)
  }
  res_SES <- df %>%
    rowwise() %>%
    mutate(
      ses_result = list(sesandpvalue_local(
        obs = .data[[obs_col]],
        rand = c_across(all_of(null_cols)),
        nreps = sum(!is.na(c_across(all_of(null_cols))))
      ))
    ) %>%
    unnest_wider(ses_result) %>%
    ungroup() %>%
    dplyr::select(usage, threat_category, Observed, SES, MeanRd, CI025Rd, CI975Rd, Pval, Nreps)
  return(res_SES)
}

################################################################################
# SHIFT IN FUNCTIONAL SPACE (FS)
################################################################################

imageTPD <- function(x, thresholdPlot = 0.99) {
  TPDList <- x$TPDc$TPDc
  imageTPD <- list()
  
  for (comm in 1:length(TPDList)) {
    percentile <- rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]] <- cbind(
      index = 1:length(TPDList[[comm]]),
      prob  = TPDList[[comm]],
      percentile
    )
    orderTPD <- order(TPDList[[comm]][, "prob"], decreasing = TRUE)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD, ]
    TPDList[[comm]][, "percentile"] <- cumsum(TPDList[[comm]][, "prob"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][, "index"]), ]
    imageTPD[[comm]] <- TPDList[[comm]]
  }
  names(imageTPD) <- names(TPDList)
  
  trait1Edges <- unique(x$data$evaluation_grid[, 1])
  trait2Edges <- unique(x$data$evaluation_grid[, 2])
  
  imageMat <- array(
    NA,
    dim = c(length(trait1Edges), length(trait2Edges), length(imageTPD)),
    dimnames = list(trait1Edges, trait2Edges, names(TPDList))
  )
  
  for (comm in 1:length(TPDList)) {
    percentileSpace <- x$data$evaluation_grid
    percentileSpace$percentile <- imageTPD[[comm]][, "percentile"]
    
    for (i in 1:length(trait2Edges)) {
      colAux <- subset(percentileSpace, percentileSpace[, 2] == trait2Edges[i])
      imageMat[, i, comm] <- colAux$percentile
    }
    
    imageMat[, , comm][imageMat[, , comm] > thresholdPlot] <- NA
  }
  
  return(imageMat)
}

plot_functional_shift_by_usage <- function(usage_name, save_path = "figures/") {
  
  message(glue::glue("Processing usage: {usage_name}"))
  traits_use  <- pca_trait$uses
  species_all <- rownames(traits_use)
  threat_vec  <- IUCN$IUCN %in% c("CR", "EN", "VU", "NT")
  used_vec    <- traits_use[[usage_name]] == 1
  comm <- matrix(
    0,
    nrow = 3,
    ncol = length(species_all),
    dimnames = list(c("ALL", "Usage", "Usagewithoutthreatened"), species_all)
  )
  comm["ALL", ] <- 1
  comm["Usage", used_vec] <- 1
  comm["Usagewithoutthreatened", used_vec & !threat_vec] <- 1
  TPDc_use <- TPDc(TPDs = TPDs_fish, sampUnit = comm)
  
  comp1 <- unique(TPDc_use$data$evaluation_grid[, 1])
  comp2 <- unique(TPDc_use$data$evaluation_grid[, 2])
  
  mat_usage     <- imageTPD(TPDc_use, thresholdPlot = 0.99)[, , "Usage"]
  mat_no_threat <- imageTPD(TPDc_use, thresholdPlot = 0.99)[, , "Usagewithoutthreatened"]
  mat_diff <- mat_usage - mat_no_threat
  
  mat_lost <- mat_usage
  mat_lost[!is.na(mat_usage) & !is.na(mat_no_threat)] <- NA
  mat_lost[!is.na(mat_usage) & is.na(mat_no_threat)]  <- 1
  
  mat_usage_full     <- imageTPD(TPDc_use, thresholdPlot = 1)[, , "Usage"]
  mat_no_threat_full <- imageTPD(TPDc_use, thresholdPlot = 1)[, , "Usagewithoutthreatened"]
  
  ncol <- 1000
  ColorRamp <- rev(scico(n = ncol, palette = "vik"))
  Min    <- -0.36
  Max    <- 0.29
  Thresh <- 0
  nHalf <- sum(!is.na(mat_diff)) / 2
  rc1 <- colorRampPalette(ColorRamp[1:500], space = "Lab")(nHalf)
  rc2 <- colorRampPalette(ColorRamp[501:1000], space = "Lab")(nHalf)
  rampcols   <- c(rc1, rc2)
  rampbreaks <- c(
    seq(Min, Thresh, length.out = nHalf + 1),
    seq(Thresh, Max, length.out = nHalf + 1)[-1]
  )
  
  cont_funspace <- contourLines(
    x = unique(TPDc_use$data$evaluation_grid[, 1]),
    y = unique(TPDc_use$data$evaluation_grid[, 2]),
    z = imageTPD(TPDc_use, thresholdPlot = 1)[, , "ALL"],
    levels = 0.999
  )
  
  cont1 <- contourLines(x = comp1, y = comp2, z = mat_usage_full,     levels = c(0.99))
  cont2 <- contourLines(x = comp1, y = comp2, z = mat_no_threat_full, levels = c(0.99))
  
  jpeg(
    filename = glue::glue("{save_path}/FS_shift_{gsub(' ', '_', usage_name)}.jpg"),
    width = 2000, height = 1600, res = 300
  )
  
  image(
    x = comp1, y = comp2, z = mat_lost,
    xlim = limX, ylim = limY,
    col = "black", breaks = c(0.5, 1.5),
    axes = FALSE, xlab = "", ylab = "", asp = 1
  )
  
  image(
    x = comp1, y = comp2, z = mat_diff,
    xlim = limX, ylim = limY,
    col = rampcols, breaks = rampbreaks,
    add = TRUE
  )
  
  for (cont in cont_funspace) {
    lines(cont$x, cont$y, lwd = 0.8, lty = 1, col = "grey30")
  }
  
  dev.off()
}

################################################################################
# DISTINCTIVENESS / MORPHOLOGICAL VARIATION
################################################################################

assign_deciles_var <- function(data, var_name = "Ui") {
  cuts <- quantile(data[[var_name]], probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  levels <- paste0("D", 1:10)
  data %>%
    mutate(Decile = cut(
      .data[[var_name]],
      breaks = cuts,
      include.lowest = TRUE,
      labels = levels
    ))
}

get_used_species_var <- function(data, var_name = "Ui") {
  data %>%
    group_by(Species) %>%
    summarise(
      !!var_name := first(.data[[var_name]]),
      Used = any(Use != "Non use" & Use_presence == 1),
      .groups = "drop"
    )
}

compute_used_proportion_var <- function(species_data, var_name = "Ui") {
  species_data %>%
    assign_deciles_var(var_name = var_name) %>%
    group_by(Decile) %>%
    summarise(
      n         = n(),
      Used_Prop = mean(Used),
      .groups   = "drop"
    )
}

bootstrap_used_proportions_var <- function(data, var_name = "Ui",
                                           n_iter = 999, prop_sample = 0.8,
                                           return_all = FALSE) {
  species_unique <- get_used_species_var(data, var_name = var_name)
  counts <- species_unique %>%
    assign_deciles_var(var_name = var_name) %>%
    count(Decile)
  res <- map_dfr(seq_len(n_iter), ~ {
    samp <- sample_frac(species_unique, prop_sample)
    compute_used_proportion_var(samp, var_name) %>%
      mutate(Iter = .x)
  })
  if (return_all) return(res)
  res %>%
    group_by(Decile) %>%
    summarise(
      n = counts$n[match(Decile, counts$Decile)],
      Mean_Used_Prop = mean(Used_Prop),
      Lower_CI       = quantile(Used_Prop, 0.025),
      Upper_CI       = quantile(Used_Prop, 0.975),
      .groups        = "drop"
    )
}

make_decile_labels_var <- function(data, var_name = "Ui") {
  cuts <- quantile(data[[var_name]], probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  labels <- paste0(
    "D", 1:10,
    " (", sprintf("%.2f", cuts[1:10]),
    "–", sprintf("%.2f", cuts[2:11]), ")"
  )
  names(labels) <- paste0("D", 1:10)
  labels
}

plot_proportions_var <- function(summary_df, original_data, var_name = "Ui") {
  labels <- make_decile_labels_var(original_data, var_name)
  overall <- original_data %>%
    get_used_species_var(var_name) %>%
    summarise(overall = mean(Used)) %>%
    pull(overall)
  
  ggplot(summary_df, aes(x = Decile, y = Mean_Used_Prop, color = Decile)) +
    geom_point(size = 4, position = position_nudge(x = 0.1)) +
    geom_errorbar(
      aes(ymin = Lower_CI, ymax = Upper_CI),
      width = 0.2,
      position = position_nudge(x = 0.1)
    ) +
    geom_hline(yintercept = overall, linetype = "dashed", color = "grey50") +
    scale_x_discrete(limits = paste0("D", 1:10), labels = labels) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    scale_color_viridis_d(begin = 0.2, end = 0.8) +
    labs(x = NULL, y = "Species used (%)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

run_statistical_tests <- function(data_long, var_name = "Ui") {
  species_unique <- get_used_species_var(data_long, var_name)
  mu <- mean(species_unique$Used)
  resamples_df <- bootstrap_used_proportions_var(data_long, var_name = var_name, return_all = TRUE)
  
  wilcox_decile_vs_mu <- resamples_df %>%
    group_by(Decile) %>%
    summarise(
      p_value = wilcox.test(Used_Prop, mu = mu)$p.value,
      median  = median(Used_Prop),
      .groups = "drop"
    )
  print(wilcox_decile_vs_mu)
  
  kr <- kruskal.test(Used_Prop ~ Decile, data = resamples_df)
  pw <- pairwise.wilcox.test(resamples_df$Used_Prop, resamples_df$Decile, p.adjust.method = "BH")
  message(sprintf("Kruskal-Wallis p = %.4f", kr$p.value))
  print(pw)
  
  summary_df <- bootstrap_used_proportions_var(data_long, var_name = var_name)
  pvals <- resamples_df %>%
    group_by(Decile) %>%
    summarise(
      p_value = wilcox.test(Used_Prop, mu = mu)$p.value,
      .groups = "drop"
    )
  
  summary_df <- summary_df %>%
    left_join(pvals, by = "Decile")
  print(summary_df)
  
  list(
    wilcox_vs_mu = wilcox_decile_vs_mu,
    kruskal      = kr,
    pairwise     = pw,
    summary      = summary_df
  )
}

compute_empirical_pvalues <- function(observed_df, null_df) {
  n_iter <- length(unique(null_df$Iter))
  observed_df %>%
    rowwise() %>%
    mutate(
      p_value = {
        q   <- Decile
        obs <- Mean_Used_Prop
        null_values <- null_df %>%
          filter(Decile == q) %>%
          pull(Used_Prop)
        (sum(null_values >= obs) + 1) / (n_iter + 1)
      }
    ) %>%
    ungroup()
}

################################################################################
# GLM MODELS (THREATENED STATUS ~ TRAIT / USE)
################################################################################

prepare_menaced_data_var <- function(data) {
  threatened_levels <- c("CR", "EN", "VU", "NT")
  
  data %>%
    filter(IUCN %in% c("CR", "EN", "VU", "NT", "LC")) %>%
    mutate(
      Menaced = if_else(IUCN %in% threatened_levels, 1, 0)
    )
}

add_use_type_var <- function(data, usage) {
  data %>%
    mutate(
      Use_type = case_when(
        .data[[usage]] == 1 ~ "Used",
        `Non use` == 1      ~ "Non use",
        TRUE                ~ NA_character_
      )
    ) %>%
    filter(!is.na(Use_type))
}

run_glm_with_nonuse_var <- function(data, usage, var_name = "Ui") {
  df <- data %>%
    add_use_type_var(usage) %>%
    prepare_menaced_data_var()
  
  if (nrow(df) < 10) return(NULL)
  
  formula_str <- as.formula(glue::glue("Menaced ~ {var_name} * Use_type"))
  
  model <- glm(formula_str, data = df, family = binomial)
  
  list(
    usage       = usage,
    variable    = var_name,
    model       = model,
    summary     = broom::tidy(model, conf.int = TRUE),
    performance = performance::model_performance(model),
    effect_plot = ggpredict(model, terms = c(var_name, "Use_type"))
  )
}

run_glm_by_usage_var <- function(data, usage, var_name = "Ui") {
  df <- data %>%
    filter(.data[[usage]] == 1) %>%
    prepare_menaced_data_var()
  
  if (nrow(df) < 10) return(NULL)
  
  formula_str <- as.formula(glue::glue("Menaced ~ {var_name}"))
  
  model <- glm(formula_str, data = df, family = binomial)
  
  list(
    usage       = usage,
    variable    = var_name,
    model       = model,
    summary     = broom::tidy(model, conf.int = TRUE),
    performance = performance::model_performance(model),
    effect_plot = ggpredict(model, terms = var_name)
  )
}

make_model_table <- function(models_list) {
  bind_rows(lapply(models_list, function(x) {
    if (is.null(x)) return(NULL)
    x$summary %>%
      filter(term == x$variable) %>%
      mutate(Usage = x$usage)
  })) %>%
    select(Usage, term, estimate, conf.low, conf.high, p.value) %>%
    arrange(p.value)
}

generate_effect_plots <- function(models_with_nonuse, variable_name = "Ui") {
  lapply(models_with_nonuse, function(x) {
    plot(x$effect_plot) +
      labs(
        title = x$usage,
        y = "Probability of being threatened",
        x = variable_name
      ) +
      theme_minimal()
  }) %>%
    patchwork::wrap_plots()
}

generate_effect_df <- function(models_with_nonuse) {
  bind_rows(lapply(models_with_nonuse, function(x) {
    df <- as.data.frame(x$effect_plot)
    df$Usage <- x$usage
    return(df)
  }))
}

get_max_val_per_usage <- function(data, usage, variable_name = "Ui") {
  data %>%
    filter(IUCN %in% c("LC", "NT", "VU", "EN", "CR")) %>%
    filter(.data[[usage]] == 1 | `Non use` == 1) %>%
    summarise(max_val = max(.data[[variable_name]], na.rm = TRUE)) %>%
    mutate(Usage = usage)
}

generate_trimmed_effect_df <- function(effects_df, raw_data, variable_name = "Ui") {
  max_vals <- bind_rows(lapply(usages, function(u) get_max_val_per_usage(raw_data, u, variable_name)))
  
  effects_df %>%
    left_join(max_vals, by = "Usage") %>%
    filter(x <= max_val)
}

plot_glm_effects <- function(effects_df_trimmed, variable_name = "Ui") {
  ggplot(effects_df_trimmed, aes(x = x, y = predicted, color = group)) +
    geom_line(size = 1.2) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high, fill = group),
      alpha = 0.25, color = NA
    ) +
    scale_color_viridis_d(option = "B", begin = 0.2, end = 0.8, name = " ") +
    scale_fill_viridis_d(option = "B", begin = 0.2, end = 0.8, name = " ") +
    facet_wrap(~Usage, scales = "free") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = glue::glue("Morphological {tolower(variable_name)}"),
      y = "Probability of being threatened",
      title = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text   = element_text(face = "bold", size = 13),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 12),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}

run_global_glm_var <- function(data, var_name = "Ui") {
  df <- data %>%
    filter(IUCN %in% c("LC", "NT", "VU", "EN", "CR")) %>%
    mutate(Menaced = if_else(IUCN %in% c("CR", "EN", "VU", "NT"), 1, 0))
  
  model <- glm(as.formula(glue::glue("Menaced ~ {var_name}")), data = df, family = binomial)
  
  list(
    model       = model,
    summary     = broom::tidy(model, conf.int = TRUE),
    performance = performance::model_performance(model),
    effect_plot = ggpredict(model, terms = var_name)
  )
}

################################################################################
# IMPUTATION ERROR (EVALUATION)
################################################################################

evaluate_imputation_phylo <- function(traitsData, traitsDataImputed, selectedTraits,
                                      meanInputed, sdInputed, traitPCA, PCAmodel,
                                      phylogeny, dimensions = 1:4, percImpute = 0.1,
                                      nboot = 100, npcoa = 5,
                                      ncores = parallel::detectCores() - 1,
                                      ref_complete_max = 1500,
                                      ntree = 30, maxiter = 2,
                                      seed = 123) {
  
  require(missForest)
  require(ape)
  require(stats)
  require(furrr)
  require(future)
  require(progressr)
  
  set.seed(seed)
  
  cat("\n[1/7] Preparation...\n")
  
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  handlers(global = TRUE)
  
  completeSp   <- rownames(na.omit(traitsData[, selectedTraits, drop = FALSE]))
  incompleteSp <- setdiff(rownames(traitsData), completeSp)
  traitWithoutNA <- traitsData[completeSp, , drop = FALSE]
  traitWithNA    <- traitsData[incompleteSp, , drop = FALSE]
  cat("   - Complete species:", nrow(traitWithoutNA), " | incomplete:", nrow(traitWithNA), "\n")
  
  traitNAmask <- traitWithNA
  if (nrow(traitNAmask) > 0) {
    traitNAmask[!is.na(traitNAmask)] <- 1
    traitNAmask[is.na(traitNAmask)]  <- NA
  }
  
  cat("[2/7] Matching phylogeny and species...\n")
  phylogeny$tip.label <- gsub("\\.", " ", phylogeny$tip.label)
  common_species <- intersect(rownames(traitsData), phylogeny$tip.label)
  phylogeny <- keep.tip(phylogeny, common_species)
  
  cat("[3/7] Phylogenetic PCoA (k =", npcoa, ")...\n")
  cophe_dist <- cophenetic.phylo(phylogeny)
  pcoa_phylo <- cmdscale(as.dist(cophe_dist), k = npcoa)
  colnames(pcoa_phylo) <- paste0("PCo", seq_len(ncol(pcoa_phylo)))
  cat("   - PCoA with", ncol(pcoa_phylo), "axes for", nrow(pcoa_phylo), "species\n")
  
  if (!all(grepl("^Comp\\.", colnames(traitPCA)))) {
    colnames(traitPCA) <- paste0("Comp.", seq_len(ncol(traitPCA)))
  }
  comp_names <- paste0("Comp.", dimensions)
  range_vals <- apply(
    traitPCA[, comp_names, drop = FALSE],
    2,
    function(x) diff(range(x, na.rm = TRUE))
  )
  
  cat("[4/7] Starting bootstrap (nboot =", nboot, ", percImpute =", percImpute, ")...\n")
  
  progressr::with_progress({
    p <- progressr::progressor(steps = nboot)
    
    results <- future_map(seq_len(nboot), function(b) {
      size <- max(1L, round(percImpute * nrow(traitWithoutNA)))
      sel_rows <- sample(rownames(traitWithoutNA), size = size, replace = FALSE)
      traitSimulNA <- traitWithoutNA[sel_rows, , drop = FALSE]
      
      if (nrow(traitNAmask) > 0) {
        randMask <- traitNAmask[
          sample(rownames(traitNAmask), nrow(traitSimulNA), replace = TRUE),
          ,
          drop = FALSE
        ]
        traitSimulNA <- traitSimulNA * randMask
      } else {
        mask <- matrix(runif(length(traitSimulNA)) < 0.1, nrow = nrow(traitSimulNA))
        traitSimulNA[mask] <- NA_real_
      }
      
      ref_pool <- setdiff(rownames(traitWithoutNA), sel_rows)
      if (length(ref_pool) > ref_complete_max) {
        ref_pool <- sample(ref_pool, ref_complete_max, replace = FALSE)
      }
      
      traitSimulAll <- rbind(
        traitSimulNA,
        traitWithNA,
        traitWithoutNA[ref_pool, , drop = FALSE]
      )
      
      pcoa_matrix <- matrix(
        NA_real_,
        nrow = nrow(traitSimulAll),
        ncol = ncol(pcoa_phylo),
        dimnames = list(rownames(traitSimulAll), colnames(pcoa_phylo))
      )
      common_pcoa_species <- intersect(rownames(traitSimulAll), rownames(pcoa_phylo))
      pcoa_matrix[common_pcoa_species, ] <- pcoa_phylo[common_pcoa_species, , drop = FALSE]
      
      traitFull_phylo <- cbind(traitSimulAll, pcoa_matrix)
      
      imputed <- missForest(
        traitFull_phylo,
        ntree = ntree,
        maxiter = maxiter,
        verbose = FALSE,
        parallelize = "no"
      )$ximp
      
      imputed_scaled <- imputed[, selectedTraits, drop = FALSE]
      m <- meanInputed[selectedTraits]
      s <- sdInputed[selectedTraits]
      imputed_scaled <- sweep(sweep(imputed_scaled, 2, m, "-"), 2, s, "/")
      
      proj     <- predict(PCAmodel, newdata = imputed_scaled[rownames(traitSimulNA), , drop = FALSE])
      proj_ref <- traitPCA[rownames(proj), comp_names, drop = FALSE]
      
      rmse <- numeric(length(dimensions))
      for (i in seq_along(dimensions)) {
        axis <- dimensions[i]
        range_axis <- range_vals[i]
        rmse[i] <- sqrt(mean((proj[, i] - proj_ref[, i])^2, na.rm = TRUE)) / range_axis
      }
      
      p(message = sprintf("iteration %d/%d done", b, nboot))
      rmse
    }, .options = furrr_options(seed = TRUE))
    
    rmse_matrix <- do.call(rbind, results)
    colnames(rmse_matrix) <- paste0("RMSE_PC", dimensions)
    
    cat("[5/7] Aggregating results...\n")
    rmse_mean <- colMeans(rmse_matrix, na.rm = TRUE) * 100
    rmse_sd   <- apply(rmse_matrix, 2, sd, na.rm = TRUE) * 100
    
    summary <- data.frame(
      Axis         = paste0("PC", dimensions),
      Mean_percent = round(as.numeric(rmse_mean), 2),
      SD_percent   = round(as.numeric(rmse_sd), 2),
      stringsAsFactors = FALSE
    )
    
    cat("[6/7] Summary (NRMSE % : mean ± SD)\n")
    print(summary)
    
    cat("[7/7] Done.\n")
    return(list(summary = summary, RMSE = rmse_matrix))
  })
}
