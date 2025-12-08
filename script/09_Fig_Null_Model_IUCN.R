# ============================================================
# Script: 09_Fig_Null_Model_IUCN.R
# Purpose: Generate figures for the IUCN-based null model of
#          functional richness loss under sequential species removal
# Author: Pierre Bouchet
# Created: 2025-12-03
# ============================================================

################################################################################
# DATA IMPORT
################################################################################

res_SES          <- readRDS("output/res_FRic_threat_SES.rds")
res_FRic_threat  <- readRDS("output/res_FRic_threat.rds")

################################################################################
# PREPARE SIMULATED & OBSERVED VALUES
################################################################################

Fric_sim <- res_FRic_threat %>%
  pivot_longer(
    cols      = starts_with("FRic_null_"),
    names_to  = "Iteration",
    values_to = "res"
  ) %>%
  rename(Usage = usage, Threat = threat_category) %>%
  dplyr::select(-FRic_obs)

Fric_obs <- res_FRic_threat %>%
  dplyr::select(usage, threat_category, FRic_obs) %>%
  distinct() %>%
  rename(Usage = usage, Threat = threat_category)

renames <- c(
  "CR"               = "-CR",
  "CR_EN"            = "-EN",
  "CR_EN_VU"         = "-VU",
  "CR_EN_VU_NT"      = "-NT",
  "CR_EN_VU_NT_DD"   = "-DD"
)

res_SES  <- res_SES  %>% mutate(threat_category = dplyr::recode(threat_category, !!!renames))
Fric_obs <- Fric_obs %>% mutate(Threat = dplyr::recode(Threat, !!!renames))
Fric_sim <- Fric_sim %>% mutate(Threat = dplyr::recode(Threat, !!!renames))

Fric_obs <- Fric_obs %>%
  left_join(
    res_SES %>% dplyr::select(Usage = usage, Threat = threat_category, Pval),
    by = c("Usage", "Threat")
  ) %>%
  mutate(point_shape = ifelse(Pval < 0.025, "p-value < 0.025", "non-significant"))

################################################################################
# CONVERT TO PERCENTAGE OF GLOBAL FRic
################################################################################

FRich_Fish  <- readRDS("output/FRich_Fish.rds")
FRic_global <- FRich_Fish["all"]

Fric_obs <- Fric_obs %>%
  mutate(res = 100 + 100 * (FRic_obs - FRich_Fish[Usage]) / FRic_global)

Fric_sim <- Fric_sim %>%
  mutate(res = 100 + 100 * (res - FRich_Fish[Usage]) / FRic_global)

################################################################################
# FILTER & ORDER
################################################################################

Fric_obs <- Fric_obs %>% filter(Usage != "Bait")
Fric_sim <- Fric_sim %>% filter(Usage != "Bait")

usage_order <- c("All uses", "Fisheries", "Aquarium", "Aquaculture", "Game fish")

Fric_obs <- Fric_obs %>% mutate(Usage = factor(Usage, levels = usage_order))
Fric_sim <- Fric_sim %>% mutate(Usage = factor(Usage, levels = usage_order))

################################################################################
# GLOBAL FIGURE
################################################################################

p <- ggplot() +
  geom_line(
    data = Fric_sim %>% arrange(Usage, Iteration, Threat),
    aes(x = Threat, y = res, group = interaction(Usage, Iteration), color = Usage),
    alpha = 0.03, size = 0.1
  ) +
  geom_point(
    data = Fric_obs,
    aes(x = Threat, y = res, color = Usage, shape = point_shape),
    size = 2.75
  ) +
  geom_line(
    data = Fric_obs,
    aes(x = Threat, y = res, group = Usage, color = Usage),
    size = 1
  ) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
  facet_wrap(~Usage, scales = "fixed", nrow = 1) +
  coord_cartesian(ylim = c(92.5, 100)) +
  scale_x_discrete(limits = c("-CR", "-EN", "-VU", "-NT", "-DD")) +
  scale_shape_manual(values = c("p-value < 0.025" = 16, "non-significant" = 1)) +
  scale_color_manual(values = c(
    "All uses"   = "#c1ea25",
    "Aquarium"   = "#63A088",
    "Game fish"  = "#D496A7",
    "Fisheries"  = "#5EB1BF",
    "Aquaculture"= "#999999"
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(
    x = " ",
    y = "Morphological richness (%)",
    color = " ",
    shape = " "
  )

print(p)

################################################################################
# PER-USAGE PANELS
################################################################################

unique_usages <- levels(Fric_obs$Usage)

for (u in unique_usages) {
  
  current_point_obs <- data.frame(
    Threat      = "Current",
    res         = 100,
    Usage       = u,
    point_shape = NA
  )
  
  obs_full <- bind_rows(
    current_point_obs,
    Fric_obs %>% filter(Usage == u)
  )
  
  current_point_sim <- Fric_sim %>%
    filter(Usage == u) %>%
    group_by(Iteration) %>%
    summarise(res = 100, Threat = "Current", .groups = "drop") %>%
    mutate(Usage = u)
  
  sim_full <- bind_rows(
    current_point_sim,
    Fric_sim %>% filter(Usage == u)
  )
  
  fig_u <- ggplot() +
    geom_line(
      data = sim_full %>% arrange(Iteration, Threat),
      aes(x = Threat, y = res, group = Iteration, color = Usage),
      alpha = 0.03, size = 0.1
    ) +
    geom_line(
      data = obs_full %>% arrange(Threat),
      aes(x = Threat, y = res, group = Usage, color = Usage),
      size = 1
    ) +
    geom_point(
      data = obs_full,
      aes(x = Threat, y = res, color = Usage, shape = point_shape),
      size = 2.75,
      na.rm = TRUE
    ) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    scale_x_discrete(limits = c("Current", "-CR", "-EN", "-VU", "-NT", "-DD")) +
    coord_cartesian(ylim = c(92.5, 100)) +
    scale_shape_manual(values = c("p-value < 0.025" = 16, "non-significant" = 1)) +
    scale_color_manual(values = c(
      "All uses"   = "#9ECF00",
      "Aquarium"   = "#63A088",
      "Game fish"  = "#D496A7",
      "Fisheries"  = "#5EB1BF",
      "Aquaculture"= "#999999"
    )) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid       = element_blank(),
      legend.position  = "none",
      strip.background = element_blank()
    ) +
    labs(
      title = u,
      x = " ",
      y = "Morphological richness (%)"
    )
  
  ggsave(
    filename = glue::glue("C:/Users/pierr/OneDrive/Documents/GitHub/fishUsages/figures/FRic_{gsub(' ', '_', u)}.jpeg"),
    plot     = fig_u,
    width    = 6,
    height   = 5,
    units    = "in",
    dpi      = 300
  )
}

################################################################################
# SAVE GLOBAL FIGURE
################################################################################

ggsave(
  filename = "plot/Fig_null_model_IUCN.png",
  plot     = p,
  dpi      = 300,
  width    = 15,
  height   = 6,
  units    = "in"
)

################################################################################
# SUMMARY TABLE IUCN × USAGE
################################################################################

IUCN <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt")
pca_trait <- readRDS("output/pca_trait.rds")
use  <- pca_trait$uses

IUCN <- IUCN %>% as.data.frame()
use  <- use  %>% as.data.frame()

if (is.null(IUCN$species)) IUCN$species <- rownames(IUCN)
if (is.null(use$species))  use$species  <- rownames(use)

merged_df <- IUCN %>%
  dplyr::select(species, IUCN) %>%
  filter(!is.na(IUCN)) %>%
  inner_join(use, by = "species")

usage_cols <- c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "Bait", "All uses")

df_long <- merged_df %>%
  pivot_longer(cols = all_of(usage_cols), names_to = "Usage", values_to = "Used") %>%
  filter(Used == 1)

summary_table <- df_long %>%
  group_by(IUCN, Usage) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

print(summary_table)
