# --- 000_functions ------------------------------------------------------------

# --- TPD ----------------------------------------------------------------------

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

# --- PCA + TPD ----------------------------------------------------------------

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

# --- Scraping -----------------------------------------------------------------

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

# --- Functional diversity -----------------------------------------------------

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

# --- Null models --------------------------------------------------------------

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

# --- SES ----------------------------------------------------------------------

sesandpvalue <- function(obs, rand, nreps, probs = c(0.025, 0.975), rnd = 3) {
  if (length(rand) < 2 || all(rand == rand[1])) {
    SES <- NA
  } else {
    SES <- (obs - mean(rand)) / sd(rand)
  }
  pValsSES <- rank(c(obs, rand), ties.method = "random")[1] / (length(rand) + 2)
  results <- round(
    c(obs, SES, mean(rand), quantile(rand, prob = probs, na.rm = TRUE), pValsSES, nreps),
    rnd
  )
  names(results) <- c("Observed", "SES", "MeanRd", "CI025Rd", "CI975Rd", "Pval", "Nreps")
  return(results)
}

get_SES <- function(obs_df, sim_df, probs = c(0.025, 0.975), rnd = 6) {
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

get_SES_from_PCA_results <- function(results_list, probs = c(0.025, 0.975), rnd = 4) {
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

# --- Functional shifts --------------------------------------------------------

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

draw_functional_shift <- function(usage_name, pca_trait, IUCN, TPDs_fish,
                                  limX = c(-7, 7), limY = c(-7, 7)) {

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
  TPDc_use <- TPD::TPDc(TPDs = TPDs_fish, sampUnit = comm)

  comp1 <- unique(TPDc_use$data$evaluation_grid[, 1])
  comp2 <- unique(TPDc_use$data$evaluation_grid[, 2])

  img_99 <- imageTPD(TPDc_use, thresholdPlot = 0.99)
  mat_usage     <- img_99[, , "Usage"]
  mat_no_threat <- img_99[, , "Usagewithoutthreatened"]
  mat_diff <- mat_usage - mat_no_threat

  mat_lost <- mat_usage
  mat_lost[!is.na(mat_usage) & !is.na(mat_no_threat)] <- NA
  mat_lost[!is.na(mat_usage) & is.na(mat_no_threat)]  <- 1

  ColorRamp <- rev(scico::scico(n = 1000, palette = "vik"))
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
    x = comp1, y = comp2,
    z = imageTPD(TPDc_use, thresholdPlot = 1)[, , "ALL"],
    levels = 0.999
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
}

draw_shift_legend <- function() {
  Min    <- -0.3
  Max    <- 0.3
  Thresh <- 0
  nHalf  <- 500

  ColorRamp <- rev(scico::scico(n = 1000, palette = "vik"))
  rc1 <- colorRampPalette(ColorRamp[1:nHalf], space = "Lab")(nHalf)
  rc2 <- colorRampPalette(ColorRamp[(nHalf + 1):1000], space = "Lab")(nHalf)
  rampbreaks <- c(
    seq(Min, Thresh, length.out = nHalf + 1),
    seq(Thresh, Max, length.out = nHalf + 1)[-1]
  )

  par(mar = c(4, 5, 2, 2))
  fields::image.plot(
    zlim = c(Min, Max), legend.only = TRUE,
    col = c(rc1, rc2), breaks = rampbreaks,
    horizontal = FALSE, legend.width = 1.2, legend.mar = 4,
    axis.args = list(at = seq(-0.3, 0.3, by = 0.1), labels = paste0(seq(-30, 30, by = 10), "%"))
  )
}

# --- Distinctiveness ----------------------------------------------------------

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

# --- Imputation error ---------------------------------------------------------

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

# --- Deficit maps -------------------------------------------------------------

imageTPD_core <- function(tpd_c, thresholdPlot = 0.999) {
  TPDList <- tpd_c$TPDc$TPDc
  percentile_tables <- vector("list", length(TPDList))
  for (k in seq_along(TPDList)) {
    tmp <- cbind(index = seq_along(TPDList[[k]]), prob = TPDList[[k]])
    tmp <- tmp[order(tmp[, "prob"], decreasing = TRUE), , drop = FALSE]
    tmp <- cbind(tmp, percentile = cumsum(tmp[, "prob"]))
    percentile_tables[[k]] <- tmp[order(tmp[, "index"]), , drop = FALSE]
  }
  xvals <- unique(tpd_c$data$evaluation_grid[, 1])
  yvals <- unique(tpd_c$data$evaluation_grid[, 2])
  out <- array(NA_real_, dim = c(length(xvals), length(yvals), length(TPDList)),
               dimnames = list(xvals, yvals, names(TPDList)))
  for (k in seq_along(TPDList)) {
    df <- tpd_c$data$evaluation_grid
    df$percentile <- percentile_tables[[k]][, "percentile"]
    for (j in seq_along(yvals)) out[, j, k] <- df[df[, 2] == yvals[j], , drop = FALSE]$percentile
    out[, , k][out[, , k] > thresholdPlot] <- NA_real_
  }
  out
}

occupancy_counts <- function(tpd_c) {
  lapply(tpd_c$TPDc$RTPDs, function(m) { m[m > 0] <- 1; as.integer(rowSums(m)) })
}

vec_to_mat <- function(v, eval_grid) {
  xvals <- unique(eval_grid[, 1]); yvals <- unique(eval_grid[, 2])
  out <- matrix(NA_real_, nrow = length(xvals), ncol = length(yvals), dimnames = list(xvals, yvals))
  tmp <- eval_grid; tmp$val <- v
  for (j in seq_along(yvals)) out[, j] <- tmp[tmp[, 2] == yvals[j], , drop = FALSE]$val
  out
}

deficit_maps <- function(TPDs, comm) {
  tpd_c <- TPD::TPDc(TPDs = TPDs, sampUnit = comm)
  list(
    core_099  = imageTPD_core(tpd_c, thresholdPlot = 0.999),
    core_full = imageTPD_core(tpd_c, thresholdPlot = 1),
    counts    = occupancy_counts(tpd_c),
    eval_grid = tpd_c$data$evaluation_grid
  )
}

draw_deficit_panel <- function(maps, catg, limX, limY, xlab, ylab, palette,
                               ncol = 1000, support_level = 0.999, contour_level = 0.999) {
  xv <- unique(maps$eval_grid[, 1]); yv <- unique(maps$eval_grid[, 2])
  core_all <- maps$core_099[, , "ALL"]

  n_all <- maps$counts[["ALL"]]
  n_use <- maps$counts[[catg]]
  deficit_vec <- 1 - ((n_all - n_use) / n_all)
  deficit_vec[!is.finite(deficit_vec)] <- NA_real_
  deficit_vec[deficit_vec < 0] <- 0
  deficit_vec[deficit_vec > 1] <- 1
  deficit_mat <- vec_to_mat(deficit_vec, maps$eval_grid)
  deficit_mat[is.na(core_all)] <- NA_real_

  ColorRamp <- palette(ncol)

  image(x = xv, y = yv, z = core_all, xlim = limX, ylim = limY, col = ColorRamp,
        xaxs = "r", yaxs = "r", axes = FALSE, asp = 1, xlab = "", ylab = "")
  for (cont in contourLines(x = xv, y = yv, z = maps$core_full[, , "ALL"], levels = support_level)) {
    polygon(x = cont$x, y = cont$y, col = 1, border = NA)
  }
  image(x = xv, y = yv, z = deficit_mat, xlim = limX, ylim = limY, col = ColorRamp,
        add = TRUE, xlab = "", ylab = "")
  for (cont in contourLines(x = xv, y = yv, z = maps$core_full[, , "ALL"], levels = contour_level)) {
    lines(cont$x, cont$y, lwd = 1.5, lty = 1, col = "black")
  }

  box(which = "plot")
  axis(1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  axis(2, las = 1, tcl = 0.3, lwd = 0.8, cex.axis = 1.1)
  mtext(xlab, side = 1, line = 2.2, cex = 1.0)
  mtext(ylab, side = 2, line = 2.6, cex = 1.0)
  title(main = catg, cex.main = 1.2)

  invisible(ColorRamp)
}

draw_deficit_row <- function(maps, usages, limX, limY, xlab, ylab, palette, ncol = 1000) {
  layout(matrix(1:6, nrow = 1), widths = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.1))
  par(mar = c(5.2, 5.2, 2.5, 0.5))
  for (catg in usages) {
    ColorRamp <- draw_deficit_panel(maps, catg, limX, limY, xlab, ylab, palette, ncol)
  }
  par(mar = c(5.2, 0.8, 2.5, 0.8))
  plot(c(0, 2), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  rasterImage(as.raster(matrix(ColorRamp, ncol = 1)), xleft = 0, ybottom = 0, xright = 1, ytop = 1)
  y_ticks <- c(1, 0.75, 0.5, 0.25, 0)
  segments(x0 = 1.00, x1 = 1.10, y0 = y_ticks, y1 = y_ticks, lwd = 1.2)
  graphics::text(x = 1.55, y = y_ticks, labels = paste0(c(100, 75, 50, 25, 0), "%"), cex = 1)
  rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, border = "black", lwd = 1)
}

# --- Trait labels -------------------------------------------------------------

trait_labels <- c(
  es  = "Relative eye size",
  ep  = "Vertical eye position",
  ms  = "Relative maxillary length",
  mp  = "Oral gape position",
  elo = "Body elongation",
  wid = "Body lateral shape",
  pp  = "Pectoral fin vertical position",
  ps  = "Pectoral fin size",
  cs  = "Caudal peduncle throttling",
  svl = "Standard body length",
  bm  = "Body mass"
)

# --- Correlation circle -------------------------------------------------------

plot_cor_circle <- function(pca_trait, pc_x = 1, pc_y = 2, stretch = 1.6) {
  loadings <- unclass(pca_trait$pca_object$loadings)
  loadings[, 1] <- -loadings[, 1]
  sdev <- pca_trait$pca_object$sdev
  pct  <- round(100 * sdev^2 / sum(sdev^2), 1)

  arrows_df <- tibble(x = stretch * loadings[, pc_x], y = stretch * loadings[, pc_y])
  circle_df <- tibble(angle = seq(0, 2 * pi, length.out = 300), x = cos(angle), y = sin(angle))

  ggplot() +
    geom_path(data = circle_df, aes(x = x, y = y), color = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +
    geom_segment(
      data = arrows_df, aes(x = 0, y = 0, xend = x, yend = y),
      color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed"), linewidth = 1.1
    ) +
    labs(
      x = paste0("PC", pc_x, " (", pct[pc_x], "%)"),
      y = paste0("PC", pc_y, " (", pct[pc_y], "%)"),
      title = " "
    ) +
    coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15)) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.title = element_text(size = 12)
    )
}

# --- Figure layout ------------------------------------------------------------

fig_font <- "Arial"

fig_save <- function(file, width, height, draw, dpi = 600) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  grDevices::cairo_pdf(paste0(file, ".pdf"), width = width / 72, height = height / 72)
  fig_draw(width, height, draw)
  grDevices::dev.off()
  ragg::agg_png(paste0(file, ".png"), width = width / 72, height = height / 72,
                units = "in", res = dpi, background = "white")
  fig_draw(width, height, draw)
  grDevices::dev.off()
  invisible(file)
}

fig_draw <- function(width, height, draw) {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::pushViewport(grid::viewport(xscale = c(0, width), yscale = c(height, 0)))
  draw()
  grid::popViewport()
}

fig_viewport <- function(x, y, width, height) {
  grid::viewport(
    x = grid::unit(x, "native"), y = grid::unit(y, "native"),
    width = grid::unit(width, "bigpts"), height = grid::unit(height, "bigpts"),
    just = c("left", "top")
  )
}

fig_plot <- function(plot, x, y, width, height) {
  vp <- fig_viewport(x, y, width, height)
  if (inherits(plot, "gg")) {
    print(plot, vp = vp, newpage = FALSE)
  } else {
    grid::pushViewport(vp)
    grid::grid.draw(plot)
    grid::popViewport()
  }
}

panel_png <- function(plot, width, height, res = 300, device = c("ragg", "png"),
                      bg = "white", ...) {
  file <- tempfile(fileext = ".png")
  if (match.arg(device) == "ragg") {
    ragg::agg_png(file, width = width, height = height, units = "px", res = res,
                  background = bg)
  } else {
    grDevices::png(file, width = width, height = height, res = res, bg = bg, ...)
  }
  if (is.function(plot)) plot() else print(plot)
  grDevices::dev.off()
  file
}

fig_image <- function(file, x, y, width, height, clip = NULL, flip_y = FALSE) {
  img <- png::readPNG(file)
  if (flip_y) img <- img[dim(img)[1]:1, , , drop = FALSE]
  if (!is.null(clip)) {
    grid::pushViewport(grid::viewport(
      x = grid::unit(clip[1], "native"), y = grid::unit(clip[2], "native"),
      width = grid::unit(clip[3] - clip[1], "bigpts"), height = grid::unit(clip[4] - clip[2], "bigpts"),
      just = c("left", "top"), clip = "on",
      xscale = c(clip[1], clip[3]), yscale = c(clip[4], clip[2])
    ))
    on.exit(grid::popViewport())
  }
  grid::grid.raster(
    img,
    x = grid::unit(x, "native"), y = grid::unit(y, "native"),
    width = grid::unit(width, "bigpts"), height = grid::unit(height, "bigpts"),
    just = c("left", "top"), interpolate = TRUE
  )
}

fig_text <- function(label, x, y, size, face = "plain", col = "black", rot = 0, hjust = 0) {
  grid::grid.text(
    label, x = grid::unit(x, "native"), y = grid::unit(y, "native"),
    hjust = hjust, vjust = 0, rot = rot,
    gp = grid::gpar(fontfamily = fig_font, fontface = face, fontsize = size, col = col)
  )
}

fig_text_runs <- function(runs, x, y, col = "black") {
  for (r in runs) {
    dy   <- if (is.null(r$dy)) 0 else r$dy
    face <- if (is.null(r$face)) "plain" else r$face
    fig_text(r[[1]], x, y + dy, r[[2]], face, col)
    gp <- grid::gpar(fontfamily = fig_font, fontface = face, fontsize = r[[2]])
    x  <- x + grid::convertWidth(grid::grobWidth(grid::textGrob(r[[1]], gp = gp)),
                                 "bigpts", valueOnly = TRUE)
  }
}

fig_texts <- function(df) {
  df <- as.data.frame(df)
  if (is.null(df$face)) df$face <- "plain"
  if (is.null(df$col))  df$col  <- "black"
  if (is.null(df$rot))  df$rot  <- 0
  if (!nrow(df)) return(invisible())
  for (i in seq_len(nrow(df))) {
    fig_text(df$label[i], df$x[i], df$y[i], df$size[i], df$face[i], df$col[i], df$rot[i])
  }
}

fig_line <- function(x, y, col = "black", lwd = 0.5, lty = "solid", fill = NA,
                     arrow = NULL, lineend = "butt") {
  gp <- grid::gpar(col = col, lwd = lwd * 96 / 72, lty = lty, fill = fill,
                   lineend = lineend, linejoin = "mitre")
  if (is.na(fill)) {
    grid::grid.lines(grid::unit(x, "native"), grid::unit(y, "native"), arrow = arrow, gp = gp)
  } else {
    grid::grid.polygon(grid::unit(x, "native"), grid::unit(y, "native"), gp = gp)
  }
}

fig_circle <- function(x, y, r, col = "black", lwd = 0.5, fill = NA) {
  grid::grid.circle(grid::unit(x, "native"), grid::unit(y, "native"), grid::unit(r, "bigpts"),
                    gp = grid::gpar(col = col, lwd = lwd * 96 / 72, fill = fill))
}

fig_arrow <- function(x0, y0, x1, y1, lwd, col = "#808080", dashed = TRUE, dot = TRUE,
                      both_ends = FALSE) {
  len <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  u   <- c(x1 - x0, y1 - y0) / len
  if (dashed) {
    s <- seq(0, len, by = 5 * lwd)
    e <- pmin(s + 4 * lwd, len)
    grid::grid.segments(
      x0 + s * u[1], y0 + s * u[2], x0 + e * u[1], y0 + e * u[2], default.units = "native",
      gp = grid::gpar(col = col, lwd = lwd * 96 / 72, lineend = "butt")
    )
  } else {
    fig_line(c(x0, x1), c(y0, y1), col = col, lwd = lwd)
  }
  head <- function(px, py, v) {
    tip <- c(px, py) + 0.25 * lwd * v
    arm <- function(a) tip - 4.24 * lwd * c(v[1] * cos(a) - v[2] * sin(a), v[1] * sin(a) + v[2] * cos(a))
    p <- rbind(arm(pi / 4), tip, arm(-pi / 4))
    fig_line(p[, 1], p[, 2], col = col, lwd = lwd)
  }
  head(x1, y1, u)
  if (both_ends) head(x0, y0, -u)
  if (dot) fig_circle(x0, y0, 2.5 * lwd, col = NA, fill = col)
}

# --- PhyloPic silhouettes -----------------------------------------------------

fetch_phylopic_cache <- function(manifest, cache_file) {
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
  new_uuid <- setdiff(unique(manifest$uuid), names(cache))
  if (length(new_uuid) > 0L) {
    for (u in new_uuid) cache[[u]] <- rphylopic::get_phylopic(u, format = "vector")
    saveRDS(cache, cache_file)
  }
  setNames(cache[manifest$uuid], manifest$species)
}

fig_silhouette <- function(img, x0, y0, x1, y1, flip = "none") {
  if (flip == "h") img <- rphylopic::flip_phylopic(img, horizontal = TRUE, vertical = FALSE)
  grid::pushViewport(fig_viewport(x0, y0, x1 - x0, y1 - y0))
  grid::grid.draw(grImport2::pictureGrob(
    img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
    xscale = img@summary@xscale, yscale = img@summary@yscale,
    distort = TRUE, expansion = 0, clip = "off"
  ))
  grid::popViewport()
}
