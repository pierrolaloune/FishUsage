# ============================================================
# Script: 000_library.R
# Purpose: Install (if needed) and load all required packages
# ============================================================

# --- PACKAGES ------------------------------------------------

required_packages <- unique(c(
  "ade4", "ape", "berryFunctions", "betapart", "biscale", "cowplot",
  "data.table", "dplyr", "funrar", "geiger", "ggplot2", "ggpubr",
  "lsmeans", "missForest", "motmot", "multcomp", "mvMORPH", "paleotree",
  "pals", "paran", "phylobase", "phytools", "picante", "plotly",
  "plotrix", "psych", "quanteda", "ratematrix", "RColorBrewer",
  "readr", "rgbif", "rnaturalearth", "rredlist", "sf", "shape",
  "stats", "tidyr", "TPD", "vegan", "VennDiagram", "viridis",
  "wesanderson", "rvest", "xml2", "stringr", "purrr", "glue",
  "furrr", "future", "progressr", "funspace", "forcats", "ggeffects",
  "patchwork", "performance", "scico", "fields", "rfishbase",
  "pbapply", "ggrepel", "AICcmodavg", "lme4", "DHARMa", "missRanger"
))

# Install missing packages
missing_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}

# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))
