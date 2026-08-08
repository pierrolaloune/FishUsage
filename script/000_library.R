# ------------------------------------------------------------------------------
# Script : 000_library
# Author : P. Bouchet
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# METHODOLOGICAL SUMMARY
# ------------------------------------------------------------------------------

# This script declares every R package used across the project, installs those
# that are missing from the local library, and loads them all.
#
# Run it first, at the start of every session, before any other script.

# ------------------------------------------------------------------------------
# Required packages
# ------------------------------------------------------------------------------

required_packages <- unique(c(
  "ade4", "ape", "berryFunctions", "betapart", "biscale", "cowplot",
  "data.table", "dplyr", "funrar", "geiger", "ggplot2", "ggpubr",
  "lsmeans", "missForest", "motmot", "multcomp", "mvMORPH", "paleotree",
  "pals", "paran", "phytools", "picante", "plotly",
  "plotrix", "psych", "quanteda", "ratematrix", "RColorBrewer",
  "readr", "rgbif", "rnaturalearth", "rredlist", "sf", "shape",
  "stats", "tidyr", "TPD", "vegan", "VennDiagram", "viridis",
  "wesanderson", "rvest", "xml2", "stringr", "purrr", "glue",
  "furrr", "future", "progressr", "funspace", "forcats", "ggeffects",
  "patchwork", "performance", "scico", "fields", "rfishbase",
  "pbapply", "ggrepel", "AICcmodavg", "lme4", "DHARMa", "missRanger",
  "paletteer", "naniar", "mice", "VIM", "visdat", "broom", "knitr", "mgcv",
  "emayili"
))

# ------------------------------------------------------------------------------
# Install and load
# ------------------------------------------------------------------------------

# ---- Install the packages that are not yet available locally ----
# The first run can take a long time; later runs skip this step entirely.
missing_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}

# ---- Load every package ----
invisible(lapply(required_packages, library, character.only = TRUE))
