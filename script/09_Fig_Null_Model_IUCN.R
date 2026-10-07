# --- 09_Fig_Null_Model_IUCN ---------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")
source("script/000_layout.R")

# --- Data ---------------------------------------------------------------------

res_SES         <- readRDS("output/res_FRic_threat_SES.rds")
res_FRic_threat <- readRDS("output/res_FRic_threat.rds")
FRich_Fish      <- readRDS("output/FRich_Fish.rds")
pca_trait       <- readRDS("output/pca_trait.rds")
TPDs_fish       <- readRDS("output/TPD_2D.rds")
uni_clean       <- readRDS("output/uni_clean.rds")
IUCN_table      <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt")

threat_labels <- c(
  "CR"             = "-CR",
  "CR_EN"          = "-EN",
  "CR_EN_VU"       = "-VU",
  "CR_EN_VU_NT"    = "-NT",
  "CR_EN_VU_NT_DD" = "-DD"
)

usage_order <- c("All uses", "Fisheries", "Aquarium", "Game fish", "Aquaculture")

use_cols <- c(
  "All uses"    = "#9ECF00",
  "Aquarium"    = "#63A088",
  "Game fish"   = "#D496A7",
  "Fisheries"   = "#5EB1BF",
  "Aquaculture" = "#999999"
)

# --- FRic trajectories --------------------------------------------------------

Fric_sim <- res_FRic_threat %>%
  pivot_longer(starts_with("FRic_null_"), names_to = "Iteration", values_to = "res") %>%
  rename(Usage = usage, Threat = threat_category) %>%
  dplyr::select(-FRic_obs) %>%
  mutate(Threat = recode(Threat, !!!threat_labels))

Fric_obs <- res_FRic_threat %>%
  dplyr::select(usage, threat_category, FRic_obs) %>%
  distinct() %>%
  rename(Usage = usage, Threat = threat_category) %>%
  mutate(Threat = recode(Threat, !!!threat_labels)) %>%
  left_join(
    res_SES %>%
      mutate(threat_category = recode(threat_category, !!!threat_labels)) %>%
      dplyr::select(Usage = usage, Threat = threat_category, Pval),
    by = c("Usage", "Threat")
  ) %>%
  mutate(point_shape = ifelse(Pval < 0.025, "p-value < 0.025", "non-significant"))

FRic_global <- FRich_Fish["all"]

Fric_obs <- Fric_obs %>% mutate(res = 100 + 100 * (FRic_obs - FRich_Fish[Usage]) / FRic_global)
Fric_sim <- Fric_sim %>% mutate(res = 100 + 100 * (res - FRich_Fish[Usage]) / FRic_global)

Fric_obs <- Fric_obs %>% filter(Usage != "Bait") %>% mutate(Usage = factor(Usage, levels = usage_order))
Fric_sim <- Fric_sim %>% filter(Usage != "Bait") %>% mutate(Usage = factor(Usage, levels = usage_order))

delta_nt <- Fric_obs %>%
  filter(Threat == "-NT") %>%
  transmute(Usage = as.character(Usage), label = sprintf(" = %.2f%%", res - 100))

# --- FRic loss ----------------------------------------------------------------

plot_fric_use <- function(u) {
  obs_full <- bind_rows(
    data.frame(Threat = "Current", res = 100, Usage = u, point_shape = NA),
    Fric_obs %>% filter(Usage == u)
  )
  sim_full <- bind_rows(
    Fric_sim %>% filter(Usage == u) %>% group_by(Iteration) %>%
      summarise(res = 100, Threat = "Current", .groups = "drop") %>% mutate(Usage = u),
    Fric_sim %>% filter(Usage == u)
  )

  ggplot() +
    geom_line(
      data = sim_full %>% arrange(Iteration, Threat),
      aes(x = Threat, y = res, group = Iteration, color = Usage),
      alpha = 0.03, linewidth = 0.1
    ) +
    geom_line(
      data = obs_full %>% arrange(Threat),
      aes(x = Threat, y = res, group = Usage, color = Usage),
      linewidth = 1
    ) +
    geom_point(
      data = obs_full,
      aes(x = Threat, y = res, color = Usage, shape = point_shape),
      size = 2.75, na.rm = TRUE
    ) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    scale_x_discrete(limits = c("Current", unname(threat_labels))) +
    coord_cartesian(ylim = c(92.5, 100)) +
    scale_shape_manual(values = c("p-value < 0.025" = 16, "non-significant" = 1)) +
    scale_color_manual(values = use_cols) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(), legend.position = "none", strip.background = element_blank()) +
    labs(title = u, x = " ", y = "Morphological richness (%)")
}

fric_png <- purrr::map(fig4_uses, function(u) {
  main <- u == "All uses"
  p <- plot_fric_use(u) +
    labs(title = " ", y = if (main) "Morphological richness (%)" else " ") +
    theme(
      axis.text.x     = element_text(colour = NA),
      axis.ticks.x    = if (main) element_line() else element_line(colour = NA),
      plot.background = element_blank()
    )
  panel_png(p, 1800, 1500, res = 300, bg = "transparent")
})

# --- Functional space shifts --------------------------------------------------

IUCN <- uni_clean %>% dplyr::select(species, IUCN)

map_png <- purrr::map(fig4_uses, ~ panel_png(
  function() draw_functional_shift(.x, pca_trait, IUCN, TPDs_fish),
  2000, 1600, res = 300, device = "png", bg = "transparent"
))

legend_png <- panel_png(draw_shift_legend, 500, 1600, res = 300, device = "png", bg = "transparent")

# --- IUCN x use ---------------------------------------------------------------

use <- as.data.frame(pca_trait$uses)
if (is.null(IUCN_table$species)) IUCN_table$species <- rownames(IUCN_table)
use$species <- rownames(use)

summary_table <- IUCN_table %>%
  dplyr::select(species, IUCN) %>%
  filter(!is.na(IUCN)) %>%
  inner_join(use, by = "species") %>%
  pivot_longer(c("Fisheries", "Aquaculture", "Aquarium", "Game fish", "All uses"),
               names_to = "Usage", values_to = "Used") %>%
  filter(Used == 1) %>%
  group_by(IUCN, Usage) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

print(summary_table)

# --- Save ---------------------------------------------------------------------

fig_save("figures/fig4", width = 511, height = 225,
         draw = function() draw_fig4(fric_png, map_png, legend_png, delta_nt))

# --- Session ------------------------------------------------------------------

sessionInfo()
