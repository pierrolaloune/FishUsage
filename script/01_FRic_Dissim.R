# --- 01_FRic_Dissim -----------------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

tpd_trait <- readRDS("output/TPDs_fish.rds")
pca_trait <- readRDS("output/pca_trait.rds")

MatriceFish <- read.csv("output/MatriceFish.csv")
colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# --- FRic [LONG] --------------------------------------------------------------

# TPDc_Fish  <- TPDc_large(TPDs = tpd_trait, sampUnit = MatriceFish)
# FRich_Fish <- Calc_FRich(TPDc = TPDc_Fish)
#
# saveRDS(TPDc_Fish,  "output/TPDc_Fish.rds")
# saveRDS(FRich_Fish, "output/FRich_Fish.rds")

TPDc_Fish  <- readRDS("output/TPDc_Fish.rds")
FRich_Fish <- readRDS("output/FRich_Fish.rds")

# --- Null model [LONG] --------------------------------------------------------

# FRic_null_results <- simulate_FRic_null(
#   n_iter          = 999,
#   original_matrix = MatriceFish,
#   TPDs_object     = tpd_trait
# ) %>%
#   dplyr::filter(Usage != "all")
#
# saveRDS(FRic_null_results, "output/FRic_null_results.rds")

FRic_null_results <- readRDS("output/FRic_null_results.rds")

# --- SES ----------------------------------------------------------------------

obs_df <- data.frame(Use = names(FRich_Fish), FRich = as.numeric(FRich_Fish)) %>%
  dplyr::filter(Use != "all")

FRic_null_SES <- get_SES(obs_df = obs_df, sim_df = FRic_null_results)
print(FRic_null_SES)

plot_SES <- plot_SES_histograms(sim_df = FRic_null_results, obs_df = obs_df)

# --- Dissimilarity [LONG] -----------------------------------------------------

# dissimilarity_result <- dissim_large(TPDc_Fish)
# saveRDS(dissimilarity_result, "output/dissimilarity_result.rds")

dissimilarity_result <- readRDS("output/dissimilarity_result.rds")

dissim_shared <- as.data.frame(as.table(as.matrix(dissimilarity_result$communities$P_shared))) %>%
  dplyr::filter(
    !Var1 %in% c("all_uses", "all", "All uses"),
    !Var2 %in% c("all_uses", "all", "All uses")
  ) %>%
  dplyr::filter(as.character(Var1) <= as.character(Var2))

plot_dissim <- ggplot(dissim_shared, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 2)), size = 4, color = "black") +
  scale_fill_gradient(low = "#FEE8C8", high = "#E34A33", name = "P_shared") +
  labs(x = "Usage", y = "Usage", title = "P_shared", subtitle = "TPD-based shared probability") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_dissim)

# --- Save ---------------------------------------------------------------------

# saveRDS(FRic_null_SES, "output/FRic_null_SES.rds")
# write.csv(FRic_null_SES, "output/FRic_null_SES.csv", row.names = FALSE)

# --- Session ------------------------------------------------------------------

sessionInfo()
