# --- 000_library --------------------------------------------------------------

# --- Packages -----------------------------------------------------------------

required_packages <- c(
  "ape",
  "phytools",
  "data.table",
  "TPD",
  "ks",
  "funspace",
  "funrar",
  "mgcv",
  "missForest",
  "vegan",
  "doParallel",
  "paran",
  "mvtnorm",
  "rvest",
  "xml2",
  "rfishbase",
  "taxize",
  "future",
  "furrr",
  "progressr",
  "scico",
  "fields",
  "rphylopic",
  "grImport2",
  "ragg",
  "png",
  "scales",
  "patchwork",
  "ggplot2",
  "readr",
  "stringr",
  "tibble",
  "purrr",
  "tidyr",
  "dplyr"
)

# --- Installation -------------------------------------------------------------

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  message("Installing: ", paste(missing_packages, collapse = ", "))
  tryCatch(
    install.packages(
      missing_packages,
      repos        = "https://cloud.r-project.org",
      type         = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"),
      dependencies = TRUE
    ),
    error = function(e) message("Installation failed: ", conditionMessage(e))
  )
}

# --- Loading ------------------------------------------------------------------

load_errors <- list()

for (pkg in required_packages) {
  tryCatch(
    suppressPackageStartupMessages(library(pkg, character.only = TRUE)),
    error = function(e) load_errors[[pkg]] <<- conditionMessage(e)
  )
}

if (length(load_errors) > 0L) {
  warning(
    "Packages not loaded (scripts using them will fail):\n",
    paste0("  - ", names(load_errors), ": ", unlist(load_errors), collapse = "\n"),
    call. = FALSE, immediate. = TRUE
  )
}
