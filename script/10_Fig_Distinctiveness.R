# ============================================================
# Script: 10_dFig_distinctiveness.R
# Purpose: Generate figures and statistical models linking
#          morphological distinctiveness to human uses
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

pca_trait   <- readRDS("output/pca_trait.rds")
uni         <- readRDS("output/uni.rds")
dist        <- readRDS("output/dist.rds")
iucn        <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt", header = TRUE)

species_uses <- pca_trait$uses %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species")

################################################################################
# PREPARE LONG DATASETS (UNIQUENESS & DISTINCTIVENESS)
################################################################################

tab <- uni %>%
  left_join(species_uses, by = "species") %>%
  left_join(iucn %>% select(species, IUCN), by = "species") %>%
  mutate(`Non use` =
           if_else(Fisheries + Aquaculture + Aquarium + `Game fish` + Bait == 0, 1, 0))

uni_long <- tab %>%
  rename(Species = species) %>%
  pivot_longer(
    c(Fisheries, Aquaculture, Aquarium, `Game fish`, Bait),
    names_to = "Use",
    values_to = "Use_presence"
  ) %>%
  mutate(
    Use         = if_else(`Non use` == 1 & Use == "Fisheries", "Non use", Use),
    Use_presence= if_else(Use == "Non use", 1, Use_presence)
  ) %>%
  distinct(Species, Use, .keep_all = TRUE)

# --- distinctiveness long table ---

tab2 <- dist %>%
  left_join(species_uses, by = "species") %>%
  left_join(iucn %>% select(species, IUCN), by = "species") %>%
  mutate(`Non use` =
           if_else(Fisheries + Aquaculture + Aquarium + `Game fish` + Bait == 0, 1, 0))

dist_long <- tab2 %>%
  rename(Species = species) %>%
  pivot_longer(
    c(Fisheries, Aquaculture, Aquarium, `Game fish`, Bait),
    names_to = "Use",
    values_to = "Use_presence"
  ) %>%
  mutate(
    Use         = if_else(`Non use` == 1 & Use == "Fisheries", "Non use", Use),
    Use_presence= if_else(Use == "Non use", 1, Use_presence)
  ) %>%
  distinct(Species, Use, .keep_all = TRUE)

################################################################################
# BOOTSTRAPPED PROPORTIONS BY VARIABLE
################################################################################

df_ui   <- bootstrap_used_proportions_var(uni_long,   var_name = "Ui")
df_dist <- bootstrap_used_proportions_var(dist_long, var_name = "Dist")

p_dist <- plot_proportions_var(df_dist, dist_long, var_name = "Dist")

dist_counts <- get_used_species_var(dist_long, "Dist") %>%
  assign_deciles_var(var_name = "Dist") %>%
  count(Decile)

################################################################################
# STATISTICAL TESTS ON DISTINCTIVENESS
################################################################################

dist_species <- get_used_species_var(dist_long, var_name = "Dist") %>%
  assign_deciles_var(var_name = "Dist")

global_p <- mean(dist_species$Used)

decile_stats <- dist_species %>%
  group_by(Decile) %>%
  summarise(
    n_total = n(),
    n_used  = sum(Used),
    prop    = n_used / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    prop_pct    = prop * 100,
    p_value_raw = map2_dbl(n_used, n_total,
                           ~ binom.test(.x, .y, p = global_p)$p.value),
    ses         = (prop - global_p) /
      sqrt(global_p * (1 - global_p) / n_total),
    p_value     = sprintf("%.3f", p_value_raw)
  ) %>%
  select(Decile, n_total, n_used, prop, prop_pct, ses, p_value)

print(decile_stats)

################################################################################
# GLOBAL MODEL: DISTINCTIVENESS → PROBABILITY OF BEING USED
################################################################################

species_dist <- get_used_species_var(dist_long, var_name = "Dist")

glm_dist_use <- glm(Used ~ Dist, data = species_dist, family = binomial)
summary(glm_dist_use)

species_dist$Used <- as.numeric(species_dist$Used)
coef_info <- summary(glm_dist_use)$coefficients
slope     <- coef_info["Dist", "Estimate"]
pval      <- coef_info["Dist", "Pr(>|z|)"]
pseudo_r2 <- pR2(glm_dist_use)["McFadden"]

label_stats <- paste0(
  "β = ", sprintf("%.3f", slope),
  " | p = ", format.pval(pval, digits = 3, eps = .001)
)

# --- plot global relationship ---

p_dist <- ggplot(species_dist, aes(x = Dist, y = Used)) +
  geom_jitter(
    aes(color = Dist),
    width = 0, height = 0.05, size = 2, alpha = 0.6
  ) +
  geom_smooth(
    method      = "glm",
    method.args = list(family = "binomial"),
    se          = TRUE,
    color       = "#2B6CB0",
    fill        = "#2B6CB0",
    alpha       = 0.25,
    linewidth   = 1.2
  ) +
  scale_color_viridis_c(option = "D", begin = 0.85, end = 0.2,
                        name = "Distinctiveness") +
  scale_y_continuous("Probability of being used",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Morphological distinctiveness") +
  annotate(
    "text", x = Inf, y = 0.1, label = label_stats,
    hjust = 1.05, vjust = -0.3, size = 4.5, color = "black"
  ) +
  theme_minimal(base_size = 12)

p_dist

################################################################################
# MODELS BY USAGE CATEGORY
################################################################################

okabe_ito <- c(
  "Fisheries"   = "#E69F00",
  "Aquaculture" = "#56B4E9",
  "Aquarium"    = "#009E73",
  "Game fish"   = "#D55E00",
  "All uses"    = "#7E1E9C"
)

make_usage_plot <- function(data, usage_name) {
  data_use <- data %>%
    filter(Use == usage_name) %>%
    mutate(Used = Use_presence)
  
  model <- glm(Used ~ Dist, data = data_use, family = binomial)
  coef_info <- summary(model)$coefficients
  slope <- coef_info["Dist", "Estimate"]
  pval  <- coef_info["Dist", "Pr(>|z|)"]
  
  label_stats <- paste0(
    "β = ", sprintf("%.3f", slope),
    " | p = ", format.pval(pval, digits = 2, eps = .001)
  )
  
  ggplot(data_use, aes(x = Dist, y = Used)) +
    geom_jitter(
      width = 0, height = 0.05, size = 1.5, alpha = 0.6,
      color = okabe_ito[[usage_name]]
    ) +
    geom_smooth(
      method      = "glm",
      method.args = list(family = "binomial"),
      se          = TRUE,
      linewidth   = 1.1,
      color       = okabe_ito[[usage_name]],
      fill        = okabe_ito[[usage_name]],
      alpha       = 0.25
    ) +
    scale_y_continuous("Probability of being used", limits = c(0, 1)) +
    scale_x_continuous("Morphological distinctiveness") +
    annotate(
      "text", x = Inf, y = 0.8, label = label_stats,
      hjust = 1.05, vjust = -0.3, size = 3, color = "black"
    ) +
    ggtitle(usage_name) +
    theme_minimal(base_size = 12)
}

make_all_use_plot <- function(data) {
  species_dist <- data %>%
    group_by(Species) %>%
    summarise(
      Dist = first(Dist),
      Used = as.numeric(any(Use != "Non use" & Use_presence == 1)),
      .groups = "drop"
    )
  
  model <- glm(Used ~ Dist, data = species_dist, family = binomial)
  coef_info <- summary(model)$coefficients
  
  slope <- coef_info["Dist", "Estimate"]
  pval  <- coef_info["Dist", "Pr(>|z|)"]
  
  label_stats <- paste0(
    "β = ", sprintf("%.3f", slope),
    " | p = ", format.pval(pval, digits = 2, eps = .001)
  )
  
  ggplot(species_dist, aes(x = Dist, y = Used)) +
    geom_jitter(
      width = 0, height = 0.05, size = 1.5, alpha = 0.6,
      color = okabe_ito[["All uses"]]
    ) +
    geom_smooth(
      method      = "glm",
      method.args = list(family = "binomial"),
      se          = TRUE,
      color       = okabe_ito[["All uses"]],
      fill        = okabe_ito[["All uses"]],
      alpha       = 0.25,
      linewidth   = 1.2
    ) +
    scale_y_continuous("Probability of being used", limits = c(0, 1)) +
    scale_x_continuous("Morphological distinctiveness") +
    annotate(
      "text", x = Inf, y = 0.1, label = label_stats,
      hjust = 1.05, vjust = -0.3, size = 4.5, color = "black"
    ) +
    ggtitle("All uses") +
    theme_minimal(base_size = 12)
}

################################################################################
# COMBINE PANELS
################################################################################

plot_glm_all_and_usages <- function(data, usage_vector) {
  p_all   <- make_all_use_plot(data)
  p_uses  <- map(usage_vector, ~make_usage_plot(data, .x))
  p_comb  <- wrap_plots(p_uses, ncol = 2)
  (p_all | p_comb) + plot_layout(widths = c(1.1, 1))
}

usages_to_plot <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish")

fig_out <- plot_glm_all_and_usages(dist_long, usages_to_plot)

################################################################################
# SAVE OUTPUTS
################################################################################

ggsave(
  filename = "figures/fig2_clean.pdf",
  plot     = fig_out,
  width    = 12,
  height   = 6.5,
  units    = "in",
  dpi      = 300
)
