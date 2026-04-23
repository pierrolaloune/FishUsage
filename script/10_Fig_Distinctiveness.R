# ------------------------------------------------------------------------------
# Script : 10_Fig_distinctiveness
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script loads PCA/use data, species-level uniqueness and distinctiveness
# metrics, and IUCN status. It reshapes uniqueness and distinctiveness into long
# formats by usage category (including a derived "Non use" group), computes
# bootstrapped proportions across metric deciles, and runs binomial models
# linking distinctiveness to probability of being used globally and by use.
# It then combines panels into Figure 2 outputs, and fits alternative GAM models
# with summary extraction for reporting.

# ------------------------------------------------------------------------------
# Data import
# ------------------------------------------------------------------------------

# ---- Inputs ----
pca_trait   <- readRDS("output/pca_trait.rds")
uni         <- readRDS("output/uni.rds")
dist        <- readRDS("output/dist.rds")
iucn        <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt", header = TRUE)

species_uses <- pca_trait$uses %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species")

# ------------------------------------------------------------------------------
# Prepare long datasets (uniqueness & distinctiveness)
# ------------------------------------------------------------------------------

tab <- uni %>%
  left_join(species_uses, by = "species") %>%
  left_join(iucn %>% select(species, IUCN), by = "species") %>%
  mutate(`Non use` =
           if_else(Fisheries + Aquaculture + Aquarium + `Game fish` == 0, 1, 0))

uni_long <- tab %>%
  rename(Species = species) %>%
  pivot_longer(
    c(Fisheries, Aquaculture, Aquarium, `Game fish`),
    names_to = "Use",
    values_to = "Use_presence"
  ) %>%
  mutate(
    Use          = if_else(`Non use` == 1 & Use == "Fisheries", "Non use", Use),
    Use_presence = if_else(Use == "Non use", 1, Use_presence)
  ) %>%
  distinct(Species, Use, .keep_all = TRUE)

# --- distinctiveness long table ---

tab2 <- dist %>%
  left_join(species_uses, by = "species") %>%
  left_join(iucn %>% select(species, IUCN), by = "species") %>%
  mutate(`Non use` =
           if_else(Fisheries + Aquaculture + Aquarium + `Game fish` == 0, 1, 0))

dist_long <- tab2 %>%
  rename(Species = species) %>%
  pivot_longer(
    c(Fisheries, Aquaculture, Aquarium, `Game fish`),
    names_to = "Use",
    values_to = "Use_presence"
  ) %>%
  mutate(
    Use          = if_else(`Non use` == 1 & Use == "Fisheries", "Non use", Use),
    Use_presence = if_else(Use == "Non use", 1, Use_presence)
  ) %>%
  distinct(Species, Use, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Bootstrapped proportions by variable
# ------------------------------------------------------------------------------

df_ui   <- bootstrap_used_proportions_var(uni_long,  var_name = "Ui")
df_dist <- bootstrap_used_proportions_var(dist_long, var_name = "Dist")

p_dist <- plot_proportions_var(df_dist, dist_long, var_name = "Dist")

dist_counts <- get_used_species_var(dist_long, "Dist") %>%
  assign_deciles_var(var_name = "Dist") %>%
  count(Decile)

# ------------------------------------------------------------------------------
# Statistical tests on distinctiveness
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Global model: distinctiveness → probability of being used
# ------------------------------------------------------------------------------

species_dist <- get_used_species_var(dist_long, var_name = "Dist")

glm_dist_use <- glm(Used ~ Dist, data = species_dist, family = binomial)
summary(glm_dist_use)

species_dist$Used <- as.numeric(species_dist$Used)
coef_info <- summary(glm_dist_use)$coefficients
slope     <- coef_info["Dist", "Estimate"]
pval      <- coef_info["Dist", "Pr(>|z|)"]

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

# ------------------------------------------------------------------------------
# Models by usage category
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------------------------

plot_glm_all_and_usages <- function(data, usage_vector) {
  p_all  <- make_all_use_plot(data)
  p_uses <- map(usage_vector, ~make_usage_plot(data, .x))
  p_comb <- wrap_plots(p_uses, ncol = 2)
  (p_all | p_comb) + plot_layout(widths = c(1.1, 1))
}

usages_to_plot <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish")

fig_out <- plot_glm_all_and_usages(dist_long, usages_to_plot)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

ggsave(
  filename = "figures/fig2_clean.pdf",
  plot     = fig_out,
  width    = 12,
  height   = 6.5,
  units    = "in",
  dpi      = 300
)

# ------------------------------------------------------------------------------
# Alternative model : GAMs
# ------------------------------------------------------------------------------

# --- label: deviance explained only ---

format_pval <- function(p, digits = 3, threshold = 0.001) {
  # p can be 0 in mgcv summaries -> display as p < threshold
  if (is.na(p)) return("p = NA")
  if (p == 0 || p < threshold) return(paste0("p < ", format(threshold, scientific = FALSE)))
  paste0("p = ", formatC(p, format = "f", digits = digits))
}

extract_gam_label <- function(gam_fit) {
  s <- summary(gam_fit)
  
  dev <- s$dev.expl
  p_smooth <- s$s.table[1, "p-value"]  # single smooth: s(Dist)
  
  paste0(
    "dev. expl. = ", sprintf("%.1f", 100 * dev), "%| ",
    format_pval(p_smooth)
  )
}

label_y_by_use <- function(usage_name) {
  if (usage_name == "All uses") return(0.80)   # top
  if (usage_name == "Fisheries") return(0.15)  # bottom
  return(0.83)                                 # unchanged for others
}

theme_no_axis_titles <- function() {
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
}

# --- Updated palette ---
custom_cols <- c(
  "All uses"    = "#A6C800",
  "Fisheries"   = "#5EB1BF",
  "Aquarium"    = "#63A088",
  "Aquaculture" = "#999999",
  "Game fish"   = "#D496A7"
)

# --- Global model: Distinctiveness -> Probability of being used (GAM) ---
species_dist <- get_used_species_var(dist_long, var_name = "Dist") %>%
  mutate(Used = as.numeric(Used))

k_default <- 8

gam_dist_use <- mgcv::gam(
  Used ~ s(Dist, k = k_default),
  data   = species_dist,
  family = binomial(link = "logit"),
  method = "REML"
)

par(mfrow = c(2, 2))
gam.check(gam_dist_use)

label_stats_global <- extract_gam_label(gam_dist_use)

p_dist <- ggplot(species_dist, aes(x = Dist, y = Used)) +
  geom_jitter(
    aes(color = Dist),
    width = 0, height = 0.05, size = 2, alpha = 0.6
  ) +
  geom_smooth(
    method      = "gam",
    formula     = y ~ s(x, k = k_default),
    method.args = list(
      family = binomial(link = "logit"),
      method = "REML"
    ),
    se        = TRUE,
    color     = custom_cols[["All uses"]],
    fill      = custom_cols[["All uses"]],
    alpha     = 0.25,
    linewidth = 1.2
  ) +
  scale_color_viridis_c(option = "D", begin = 0.85, end = 0.2,
                        name = "Distinctiveness") +
  scale_y_continuous("Probability of being used",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous("Morphological distinctiveness") +
  annotate(
    "text", x = Inf, y = 0.1, label = label_stats_global,
    hjust = 1.05, vjust = -0.3, size = 4.2, color = "black"
  ) +
  theme_minimal(base_size = 12)

p_dist

# --- Models by usage category (GAM) ---
make_usage_plot_gam <- function(data, usage_name, k = 10) {
  data_use <- data %>%
    filter(Use == usage_name) %>%
    mutate(Used = Use_presence)
  
  gam_fit <- mgcv::gam(
    Used ~ s(Dist, k = k),
    data   = data_use,
    family = binomial(link = "logit"),
    method = "REML"
  )
  
  label_stats <- extract_gam_label(gam_fit)
  y_lab <- label_y_by_use(usage_name)
  
  ggplot(data_use, aes(x = Dist, y = Used)) +
    geom_jitter(
      width = 0, height = 0.05, size = 1.5, alpha = 0.6,
      color = custom_cols[[usage_name]]
    ) +
    geom_smooth(
      method      = "gam",
      formula     = y ~ s(x, k = k),
      method.args = list(
        family = binomial(link = "logit"),
        method = "REML"
      ),
      se        = TRUE,
      linewidth = 1.1,
      color     = custom_cols[[usage_name]],
      fill      = custom_cols[[usage_name]],
      alpha     = 0.25
    ) +
    scale_y_continuous(NULL, limits = c(0, 1)) +
    scale_x_continuous(NULL) +
    annotate(
      "text", x = Inf, y = y_lab, label = label_stats,
      hjust = 1.05, vjust = 1.1, size = 3, color = "black"
    ) +
    ggtitle(usage_name) +
    theme_minimal(base_size = 12)
}

make_all_use_plot_gam <- function(data, k = 10) {
  species_dist_all <- data %>%
    group_by(Species) %>%
    summarise(
      Dist = first(Dist),
      Used = as.numeric(any(Use != "Non use" & Use_presence == 1)),
      .groups = "drop"
    )
  
  gam_fit <- mgcv::gam(
    Used ~ s(Dist, k = k),
    data   = species_dist_all,
    family = binomial(link = "logit"),
    method = "REML"
  )
  
  label_stats <- extract_gam_label(gam_fit)
  
  ggplot(species_dist_all, aes(x = Dist, y = Used)) +
    geom_jitter(
      width = 0, height = 0.05, size = 1.5, alpha = 0.6,
      color = custom_cols[["All uses"]]
    ) +
    geom_smooth(
      method      = "gam",
      formula     = y ~ s(x, k = k),
      method.args = list(
        family = binomial(link = "logit"),
        method = "REML"
      ),
      se        = TRUE,
      color     = custom_cols[["All uses"]],
      fill      = custom_cols[["All uses"]],
      alpha     = 0.25,
      linewidth = 1.2
    ) +
    scale_y_continuous("Probability of being used", limits = c(0, 1)) +
    scale_x_continuous("Dist") +
    annotate(
      "text", x = Inf, y = 0.95, label = label_stats,
      hjust = 1.05, vjust = 1.1, size = 4.2, color = "black"
    ) +
    ggtitle("All uses") +
    theme_minimal(base_size = 12)
}

# --- Combine panels (UPDATED ORDER) ---
plot_gam_all_and_usages <- function(data, usage_vector, k = 10) {
  p_all  <- make_all_use_plot_gam(data, k = k)
  p_uses <- map(usage_vector, ~ make_usage_plot_gam(data, .x, k = k))
  p_comb <- wrap_plots(p_uses, ncol = 2)
  (p_all | p_comb) + plot_layout(widths = c(1.1, 1))
}

# ORDER requested:
usages_to_plot <- c("Fisheries", "Aquarium", "Game fish", "Aquaculture")

fig_out <- plot_gam_all_and_usages(dist_long, usages_to_plot, k = k_default)
fig_out

out_pdf <- "figures/Clean/fig2bis_clean.pdf"
out_png <- "figures/Clean/fig2bis_clean.png"

ggsave(
  filename = out_pdf,
  plot     = fig_out,
  device   = cairo_pdf,     # high-quality PDF with embedded fonts
  width    = 11.69,         # A4 landscape width in inches
  height   = 8.27,          # A4 landscape height in inches
  units    = "in",
  dpi      = 300,           # ignored for vector PDF, but harmless
  bg       = "white"
)
dev.off()

ggsave(
  filename = out_png,
  plot     = fig_out,
  width    = 11.69,
  height    = 8.27,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# --- GAM summary table (All uses + each usage, excluding "Non use") ---

# ---- Helpers ----
format_pval <- function(p, digits = 3, threshold = 0.001) {
  # mgcv can return p = 0 in summary tables; display nicely
  if (is.na(p)) return(NA_character_)
  if (p == 0 || p < threshold) return(paste0("<", format(threshold, scientific = FALSE)))
  formatC(p, format = "f", digits = digits)
}

make_species_level_all_uses <- function(data) {
  data %>%
    group_by(Species) %>%
    summarise(
      Dist = first(Dist),
      Used = as.numeric(any(Use != "Non use" & Use_presence == 1)),
      .groups = "drop"
    )
}

make_species_level_by_usage <- function(data, usage_name) {
  data %>%
    filter(Use == usage_name) %>%
    group_by(Species) %>%
    summarise(
      Dist = first(Dist),
      Used = as.numeric(any(Use_presence == 1)),
      .groups = "drop"
    )
}

fit_gam_extract <- function(df_species, usage_name, k = 8) {
  gam_fit <- mgcv::gam(
    Used ~ s(Dist, k = k),
    data   = df_species,
    family = binomial(link = "logit"),
    method = "REML"
  )
  
  s <- summary(gam_fit)
  
  edf   <- unname(s$s.table[1, "edf"])
  chisq <- unname(s$s.table[1, "Chi.sq"])
  pval  <- unname(s$s.table[1, "p-value"])
  dev   <- unname(s$dev.expl)
  
  tibble(
    Usage    = usage_name,
    edf      = edf,
    Chi.sq   = chisq,
    p_value  = format_pval(pval),
    dev_expl = 100 * dev
  )
}

# ---- Main: build table ----
k_default <- 8

usages_to_summarise <- sort(setdiff(unique(dist_long$Use), c("Non use", "Bait")))

df_all <- make_species_level_all_uses(dist_long)

tab_gam <- bind_rows(
  fit_gam_extract(df_all, usage_name = "All uses", k = k_default),
  map_dfr(usages_to_summarise, ~{
    df_u <- make_species_level_by_usage(dist_long, .x)
    fit_gam_extract(df_u, usage_name = .x, k = k_default)
  })
) %>%
  mutate(
    dev_expl = round(dev_expl, 1),
    edf      = round(edf, 3),
    Chi.sq   = round(Chi.sq, 1)
  ) %>%
  arrange(match(Usage, c("All uses", usages_to_summarise)))

tab_gam
