# --- 00_make_figures ----------------------------------------------------------

{
  if (!dir.exists("script")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

  figures <- c(
    "script/08_Fig_FSpaces_Usages.R",
    "script/12_FS_Shifts_TPD.R",
    "script/10_Fig_Distinctiveness.R",
    "script/09_Fig_Null_Model_IUCN.R",
    "script/11_Fig_Loadings_PCA.R",
    "script/05_Distinctiveness_IUCN.R"
  )

  for (script in figures) {
    message("Running ", script)
    source(script, local = new.env(), encoding = "UTF-8")
  }
}
