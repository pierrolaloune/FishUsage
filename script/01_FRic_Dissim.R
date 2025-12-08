# ============================================================
# Script: 01_Fric_Dissim.R
# Purpose: Compute functional richness (FRic), null models,
#          standardized effect sizes (SES), and TPD-based
#          dissimilarity among human uses
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

tpd_trait   <- readRDS("output/TPDs_fish.rds")
pca_trait   <- readRDS("output/PCA_fish.rds")
MatriceFish <- read.csv("output/MatriceFish.csv")

colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

################################################################################
# COMMUNITY TPD AND FRIC
################################################################################

TPDc_Fish <- TPDc_large(TPDs = tpd_trait, sampUnit = MatriceFish)
FRich_Fish <- Calc_FRich(TPDc = TPDc_Fish)

# saveRDS(TPDc_Fish, "output/TPDc_Fish.rds")
# saveRDS(FRich_Fish, "output/FRich_Fish.rds")

TPDc_Fish   <- readRDS("output/TPDc_Fish.rds")
FRich_Fish  <- readRDS("output/FRich_Fish.rds")

################################################################################
# NULL MODEL FOR FRIC
################################################################################

FRic_null_results <- simulate_FRic_null( # long recommended skip
  n_iter          = 999,
  original_matrix = MatriceFish,
  TPDs_object     = tpd_trait
) %>%
  dplyr::filter(Usage != "all")

# saveRDS(FRic_null_results, "output/FRic_null_results.rds")
FRic_null_results <- readRDS("output/FRic_null_results.rds")

################################################################################
# SES CALCULATION
################################################################################

obs_df <- data.frame(
  Use   = names(FRich_Fish),
  FRich = as.numeric(FRich_Fish),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(Use != "all")

FRic_null_SES <- get_SES(obs_df = obs_df, sim_df = FRic_null_results) # Table S1a

# saveRDS(FRic_null_SES, "output/FRic_null_SES.rds")
# write.csv(FRic_null_SES, "output/FRic_null_SES.csv", row.names = FALSE)

################################################################################
# SES PLOTS
################################################################################

plot_SES <- plot_SES_histograms( 
  sim_df = FRic_null_results,
  obs_df = obs_df
)

# ggsave("figures/Histogram_FRic_SES.png", plot_SES, width = 10, height = 6, dpi = 300)

################################################################################
# DISSIMILARITY AMONG USES
################################################################################

dissimilarity_result <- dissim_large(TPDc_Fish)
# saveRDS(dissimilarity_result, "output/dissimilarity_result.rds")
dissimilarity_result <- readRDS("output/dissimilarity_result.rds")

dissim_df_shared <- as.data.frame(as.table(
  as.matrix(dissimilarity_result$communities$P_shared)
))

filtered_dissim_df_shared <- dissim_df_shared %>% # Table 1
  dplyr::filter(
    !Var1 %in% c("all_uses", "all", "All uses"),
    !Var2 %in% c("all_uses", "all", "All uses")
  ) %>%
  dplyr::mutate(
    Var1 = dplyr::case_when(
      Var1 == "UsedasBait"          ~ "Bait",
      Var1 == "GameFish"            ~ "Recreational",
      Var1 == "UsedforAquaculture"  ~ "Aquaculture",
      Var1 == "Importance"          ~ "Fisheries",
      TRUE ~ Var1
    ),
    Var2 = dplyr::case_when(
      Var2 == "UsedasBait"          ~ "Bait",
      Var2 == "GameFish"            ~ "Recreational",
      Var2 == "UsedforAquaculture"  ~ "Aquaculture",
      Var2 == "Importance"          ~ "Fisheries",
      TRUE ~ Var2
    )
  ) %>%
  dplyr::filter(Var1 != "All uses", Var2 != "All uses") %>%
  dplyr::filter(as.character(Var1) <= as.character(Var2))

plot_dissim <- ggplot2::ggplot(
  filtered_dissim_df_shared,
  ggplot2::aes(x = Var1, y = Var2, fill = Freq)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = round(Freq, 2)), size = 4, color = "black") +
  ggplot2::scale_fill_gradient(
    low  = "#FEE8C8",
    high = "#E34A33",
    name = "P_shared"
  ) +
  ggplot2::labs(
    x = "Usage",
    y = "Usage",
    title = "P_shared",
    subtitle = "TPD-based shared probability"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# ggsave("figures/plot_dissim.png", plot_dissim, width = 10, height = 6, dpi = 300)
