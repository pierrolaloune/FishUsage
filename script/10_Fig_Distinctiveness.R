# --- 10_Fig_Distinctiveness ---------------------------------------------------

source("script/000_library.R")
source("script/000_functions.R")
source("script/000_layout.R")

# --- Data ---------------------------------------------------------------------

pca_trait <- readRDS("output/pca_trait.rds")
dist      <- readRDS("output/dist.rds")
iucn      <- read.table("dataPrepared/Fish/traitsWithPCOAIUCN.txt", header = TRUE)

species_uses <- pca_trait$uses %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species")

USAGE_COLS     <- c("Fisheries", "Aquarium", "Aquaculture", "Game fish")
usages_to_plot <- c("Fisheries", "Aquarium", "Aquaculture", "Game fish")

custom_cols <- c(
  "All uses"    = "#A6C800",
  "Fisheries"   = "#5EB1BF",
  "Aquarium"    = "#63A088",
  "Aquaculture" = "#999999",
  "Game fish"   = "#D496A7"
)

k_default <- 8

# --- Species tables -----------------------------------------------------------

tab2 <- dist %>%
  left_join(species_uses, by = "species") %>%
  left_join(iucn %>% dplyr::select(species, IUCN), by = "species") %>%
  mutate(`Non use` = if_else(Fisheries + Aquaculture + Aquarium + `Game fish` == 0, 1, 0))

dist_long <- tab2 %>%
  rename(Species = species) %>%
  pivot_longer(c(Fisheries, Aquaculture, Aquarium, `Game fish`), names_to = "Use", values_to = "Use_presence") %>%
  mutate(
    Use          = if_else(`Non use` == 1 & Use == "Fisheries", "Non use", Use),
    Use_presence = if_else(Use == "Non use", 1, Use_presence)
  ) %>%
  distinct(Species, Use, .keep_all = TRUE)

species_wide <- tab2 %>%
  rename(Species = species) %>%
  mutate(`All uses` = as.numeric(Fisheries + Aquaculture + Aquarium + `Game fish` > 0))

make_model_data <- function(data, usage_name) {
  data %>% transmute(Species = Species, Dist = Dist, Used = as.numeric(.data[[usage_name]] == 1))
}

check_pool <- map(c("All uses", USAGE_COLS), ~ {
  d <- make_model_data(species_wide, .x)
  tibble(Usage = .x, n = nrow(d), n_used = sum(d$Used))
}) %>%
  list_rbind()

print(check_pool)

# --- Deciles ------------------------------------------------------------------

df_dist     <- bootstrap_used_proportions_var(dist_long, var_name = "Dist")
p_prop_dist <- plot_proportions_var(df_dist, dist_long, var_name = "Dist")

dist_species <- get_used_species_var(dist_long, var_name = "Dist") %>%
  assign_deciles_var(var_name = "Dist")

global_p <- mean(dist_species$Used)

decile_stats <- dist_species %>%
  group_by(Decile) %>%
  summarise(n_total = n(), n_used = sum(Used), prop = n_used / n_total, .groups = "drop") %>%
  mutate(
    prop_pct    = prop * 100,
    p_value_raw = map2_dbl(n_used, n_total, ~ binom.test(.x, .y, p = global_p)$p.value),
    ses         = (prop - global_p) / sqrt(global_p * (1 - global_p) / n_total),
    p_value     = sprintf("%.3f", p_value_raw)
  ) %>%
  select(Decile, n_total, n_used, prop, prop_pct, ses, p_value)

print(decile_stats)

# --- GLM ----------------------------------------------------------------------

glm_dist_use <- glm(Used ~ Dist, data = make_model_data(species_wide, "All uses"), family = binomial)
summary(glm_dist_use)

# --- GAM ----------------------------------------------------------------------

gam_dist_use <- mgcv::gam(
  Used ~ s(Dist, k = k_default),
  data   = make_model_data(species_wide, "All uses"),
  family = binomial(link = "logit"),
  method = "REML"
)

par(mfrow = c(2, 2))
mgcv::gam.check(gam_dist_use)
par(mfrow = c(1, 1))

fit_gam_extract <- function(data, usage_name, k = 8) {
  df_species <- make_model_data(data, usage_name)
  gam_fit <- mgcv::gam(Used ~ s(Dist, k = k), data = df_species, family = binomial(link = "logit"), method = "REML")
  s <- summary(gam_fit)
  tibble(
    Usage    = usage_name,
    n        = s$n,
    n_used   = sum(df_species$Used),
    edf      = unname(s$s.table[1, "edf"]),
    ref_df   = unname(s$s.table[1, "Ref.df"]),
    Chi.sq   = unname(s$s.table[1, "Chi.sq"]),
    p_value  = unname(s$s.table[1, "p-value"]),
    dev_expl = 100 * unname(s$dev.expl)
  )
}

tab_gam_raw <- bind_rows(
  fit_gam_extract(species_wide, "All uses", k = k_default),
  map_dfr(USAGE_COLS, ~ fit_gam_extract(species_wide, .x, k = k_default))
)

print(tab_gam_raw, n = Inf, width = Inf)

tab_gam <- tab_gam_raw %>%
  mutate(edf = round(edf, 3), ref_df = round(ref_df, 3), Chi.sq = round(Chi.sq, 1), dev_expl = round(dev_expl, 2))

print(tab_gam)

# --- Figure 3 -----------------------------------------------------------------

format_pval_label <- function(p, digits = 3, threshold = 0.001) {
  if (is.na(p)) return("P = NA")
  if (p == 0 || p < threshold) return(paste0("P < ", format(threshold, scientific = FALSE)))
  paste0("P = ", formatC(p, format = "f", digits = digits))
}

extract_gam_label <- function(gam_fit, p_label = NA) {
  s <- summary(gam_fit)
  if (is.na(p_label)) p_label <- format_pval_label(s$s.table[1, "p-value"])
  paste0("dev. expl. = ", sprintf("%.1f", 100 * s$dev.expl), "%; ", p_label)
}

p_label_override <- c("Fisheries" = "P = 0.005")

fig_width  <- 250 / 72
fig_height <- 177 / 72
fig_scale  <- fig_width / 11.69

fig_title_size <- 5
fig_label_size <- c(main = 5, use = 3.7)

jitter_seed <- 42

mm_to_linewidth <- function(mm) mm / 25.4 * 96 / .pt

theme_fig3 <- function() {
  half_line <- 12 * fig_scale / 2
  theme_minimal(base_size = 12 * fig_scale, base_family = fig_font) +
    theme(
      plot.margin     = margin(0.58, half_line, half_line, half_line),
      plot.title      = element_text(size = fig_title_size, face = "bold", hjust = 0, margin = margin(b = 2.3, l = -7)),
      axis.title      = element_text(size = fig_title_size),
      axis.line       = element_line(colour = "black", linewidth = mm_to_linewidth(0.2)),
      plot.background = element_blank()
    )
}

theme_fig3_page <- function() {
  m <- 5.5 * fig_scale
  theme(plot.margin = margin(m, m, m, m), plot.background = element_rect(fill = "white", colour = NA))
}

gam_panel <- function(data_use, usage_name, k, linewidth) {
  ggplot(data_use, aes(x = Dist, y = Used)) +
    geom_point(
      position = position_jitter(width = 0, height = 0.05, seed = jitter_seed),
      size = 1.5 * fig_scale, stroke = 0.5 * fig_scale, alpha = 0.6,
      color = custom_cols[[usage_name]]
    ) +
    geom_smooth(
      method      = "gam",
      formula     = y ~ s(x, k = k),
      method.args = list(family = binomial(link = "logit"), method = "REML"),
      se          = TRUE,
      color       = custom_cols[[usage_name]],
      fill        = custom_cols[[usage_name]],
      alpha       = 0.25,
      linewidth   = linewidth * fig_scale
    )
}

gam_label <- function(data_use, usage_name, k, size) {
  gam_fit   <- mgcv::gam(Used ~ s(Dist, k = k), data = data_use, family = binomial(link = "logit"), method = "REML")
  label_pos <- fig3_label_pos[[usage_name]]
  annotate(
    "text", x = I(label_pos[["x"]]), y = I(label_pos[["y"]]),
    label = extract_gam_label(gam_fit, unname(p_label_override[usage_name])),
    hjust = 0, vjust = 0, size = size / .pt, family = fig_font, color = "black"
  )
}

make_usage_plot_gam <- function(data, usage_name, tag, k = 8) {
  data_use <- make_model_data(data, usage_name)
  gam_panel(data_use, usage_name, k, linewidth = 1.1) +
    scale_y_continuous(NULL, limits = c(0, 1)) +
    scale_x_continuous(NULL) +
    gam_label(data_use, usage_name, k, fig_label_size[["use"]]) +
    coord_cartesian(clip = "off") +
    ggtitle(paste(tag, usage_name)) +
    theme_fig3()
}

make_all_use_plot_gam <- function(data, tag = "a", k = 8) {
  data_all <- make_model_data(data, "All uses")
  gam_panel(data_all, "All uses", k, linewidth = 1.2) +
    scale_y_continuous("Probability of being used", limits = c(0, 1)) +
    scale_x_continuous("Morphological distinctiveness") +
    gam_label(data_all, "All uses", k, fig_label_size[["main"]]) +
    coord_cartesian(clip = "off") +
    ggtitle(paste(tag, "All uses")) +
    theme_fig3() +
    theme(
      plot.margin  = margin(0.58, 12 * fig_scale / 2, 12 * fig_scale / 2, 0),
      axis.title.x = element_text(hjust = 0.656, margin = margin(t = 1.11, b = -1.68)),
      axis.title.y = element_text(hjust = 0.505, margin = margin(r = 1.57))
    )
}

p_all  <- make_all_use_plot_gam(species_wide, tag = "a", k = k_default)
p_uses <- map2(usages_to_plot, letters[seq_along(usages_to_plot) + 1],
               ~ make_usage_plot_gam(species_wide, .x, tag = .y, k = k_default))

fig_out <- (p_all | wrap_plots(p_uses, ncol = 2)) +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(theme = theme_fig3_page())

print(fig_out)

# --- Save ---------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

ggsave("figures/fig3.pdf", fig_out, device = cairo_pdf, width = fig_width, height = fig_height, units = "in")
ggsave("figures/fig3.png", fig_out, width = fig_width, height = fig_height, units = "in", dpi = 600, bg = "white")

# --- Session ------------------------------------------------------------------

sessionInfo()
