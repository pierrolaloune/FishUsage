# --- 00_RunAll ----------------------------------------------------------------

{
  if (!dir.exists("script")) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

  pipeline <- c(
    "script/000_ScrappingData.R",
    "script/000_LoadDataR.R",
    "script/01_FRic_Dissim.R",
    "script/02_FSpaces_Usages.R",
    "script/03_PCA_mean_trait_value.R",
    "script/04_Null_model_IUCN.R",
    "script/05_Distinctiveness_IUCN.R",
    "script/06_Shift_FS.R",
    "script/07_imputation_error.R",
    "script/08_Fig_FSpaces_Usages.R",
    "script/09_Fig_Null_Model_IUCN.R",
    "script/10_Fig_Distinctiveness.R",
    "script/11_Fig_Loadings_PCA.R",
    "script/12_FS_Shifts_TPD.R",
    "script/13_Single_vs_MI.R",
    "script/100_Imputation_SI_MI.R",
    "script/web scrapping percentage.R"
  )

  for (script in pipeline) {
    message("Running ", script)
    source(script, local = new.env(), encoding = "UTF-8")
  }
}
