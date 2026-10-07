# --- 03_PCA_mean_trait_value --------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/PCA_fish.rds")

MatriceFish <- read.csv("output/MatriceFish.csv")
colnames(MatriceFish)[-1] <- gsub("\\.", " ", colnames(MatriceFish)[-1])
rownames(MatriceFish)     <- MatriceFish$X
MatriceFish$X             <- NULL

# --- Null model [LONG] --------------------------------------------------------

# resultats_null <- generate_null_means(
#   pca_trait      = pca_trait,
#   MatriceFish    = MatriceFish,
#   nb_simulations = 999
# )
# resultats_null$all <- NULL
#
# saveRDS(resultats_null, "output/PCA_mean_trait_values_results.rds")

PCA_mean_trait_values_results <- readRDS("output/PCA_mean_trait_values_results.rds")

for (usage in names(PCA_mean_trait_values_results)) {
  PCA_mean_trait_values_results[[usage]]$observed["Comp.1"]    <- -PCA_mean_trait_values_results[[usage]]$observed["Comp.1"]
  PCA_mean_trait_values_results[[usage]]$simulated[, "Comp.1"] <- -PCA_mean_trait_values_results[[usage]]$simulated[, "Comp.1"]
}

# --- SES ----------------------------------------------------------------------

PCA_mean_trait_values_SES <- get_SES_from_PCA_results(results_list = PCA_mean_trait_values_results)
print(PCA_mean_trait_values_SES)

# --- Null distributions -------------------------------------------------------

df_plot <- purrr::map_dfr(names(PCA_mean_trait_values_results), function(usage) {
  observed <- PCA_mean_trait_values_results[[usage]]$observed[c("Comp.1", "Comp.2")]
  as.data.frame(PCA_mean_trait_values_results[[usage]]$simulated)[, c("Comp.1", "Comp.2")] %>%
    pivot_longer(everything(), names_to = "Component", values_to = "Simulated_value") %>%
    mutate(Usage = usage, Observed_value = observed[Component])
}) %>%
  mutate(Usage = factor(Usage, levels = c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "Bait", "All uses")))

p <- ggplot(df_plot, aes(x = Simulated_value)) +
  geom_histogram(bins = 50, fill = "#69b3a2", color = "black") +
  geom_vline(aes(xintercept = Observed_value), color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_grid(Usage ~ Component, scales = "free") +
  labs(
    x = "Simulated PCA score", y = "Frequency",
    title = "Simulated PCA distributions vs Observed values",
    subtitle = "PCA Components 1 and 2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text    = element_text(face = "bold"),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p)

# --- Save ---------------------------------------------------------------------

# saveRDS(PCA_mean_trait_values_SES, "output/PCA_mean_trait_values_SES.rds")

# --- Session ------------------------------------------------------------------

sessionInfo()
