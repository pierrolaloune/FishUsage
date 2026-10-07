# FishUsage

[![Paper](https://img.shields.io/badge/Nature%20Communications-10.1038%2Fs41467--026--78568--9-blue)](https://doi.org/10.1038/s41467-026-78568-9)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21873314.svg)](https://doi.org/10.5281/zenodo.21873314)

Data and code for the paper:

> Bouchet, P., Brosse, S. & Toussaint, A. **Human targeting of
> morphologically unique fishes amplifies the risk of functional erosion.**
> *Nature Communications* (2026).
> [https://doi.org/10.1038/s41467-026-78568-9](https://doi.org/10.1038/s41467-026-78568-9)

The analysis covers **8,970 freshwater fish species** and five categories of human
use (fisheries, aquaculture, aquarium, bait, game fish). It asks where used
species sit in the morphological space of freshwater fishes, whether
morphologically distinctive species are more likely to be targeted, and how much
of that space would be lost if threatened species disappeared.

## Repository contents

```
FishUsage/
├── script/                            # Analysis pipeline
│   ├── 00_RunAll.R                    # Runs the whole pipeline
│   ├── 00_make_figures.R              # Draws the figures of the paper
│   ├── 000_library.R                  # Packages
│   ├── 000_functions.R                # Custom functions
│   ├── 000_layout.R                   # Page layout of the figures
│   ├── 000_LoadDataR.R                # Traits, phylogeny, IUCN status, human uses
│   ├── 000_ScrappingData.R            # Human uses scraped from FishBase pages
│   ├── 01_FRic_Dissim.R               # Functional richness and dissimilarity
│   ├── 02_FSpaces_Usages.R            # Functional spaces per use category
│   ├── 03_PCA_mean_trait_value.R      # Null model on mean PCA scores
│   ├── 04_Null_model_IUCN.R           # FRic loss under nested threat scenarios
│   ├── 05_Distinctiveness_IUCN.R      # Uniqueness and distinctiveness (Fig. S5)
│   ├── 06_Shift_FS.R                  # 2D TPD used by the shift maps (Fig. 4)
│   ├── 07_imputation_error.R          # missForest imputation error
│   ├── 08_Fig_FSpaces_Usages.R        # Functional spaces (Figs 1, S1)
│   ├── 09_Fig_Null_Model_IUCN.R       # Functional richness loss (Fig. 4)
│   ├── 10_Fig_Distinctiveness.R       # Distinctiveness (Fig. 3)
│   ├── 11_Fig_Loadings_PCA.R          # PCA loadings (Fig. S2)
│   ├── 12_FS_Shifts_TPD.R             # Functional deficit maps (Figs 2, S3)
│   ├── 13_Single_vs_MI.R              # Single vs multiple imputation
│   ├── 100_Imputation_SI_MI.R         # 100 missForest imputations
│   └── web scrapping percentage.R     # Contribution of the scraping over rfishbase
├── dataPrepared/                      # Prepared input data
├── output/                            # Precomputed results (~190 MB)
├── figures/                           # Figures 1-4 and S1-S5 (PDF, PNG 600 dpi)
└── FishUsage.Rproj                    # RStudio project
```

## Data

All data needed to reproduce the analyses and figures are included in
`dataPrepared/` (traits, phylogeny, IUCN categories, human uses, phylogenetic
PCoA) and `output/` (precomputed results). They are also archived on Zenodo
([10.5281/zenodo.21873314](https://doi.org/10.5281/zenodo.21873314)).

| Dataset              | Source                                                                                                                                                                                                                  |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Morphological traits | Brosse, S. et al. FISHMORPH: a global database on morphological traits of freshwater fishes. *Glob. Ecol. Biogeogr.* **30**, 2330–2336 (2021). [doi:10.6084/m9.figshare.14891412](https://doi.org/10.6084/m9.figshare.14891412) |
| Phylogeny            | Rabosky, D. L. et al. Data from: An inverse latitudinal gradient in speciation rate for marine fishes. Dryad (2019). [doi:10.5061/dryad.fc71cp4](https://doi.org/10.5061/dryad.fc71cp4)                                  |
| Human uses           | Froese, R. & Pauly, D. FishBase (2025), via [`rfishbase`](https://docs.ropensci.org/rfishbase/) and the species summary pages of [fishbase.org](https://www.fishbase.org/)                                                |
| Conservation status  | IUCN Red List of Threatened Species, 2024 ([iucnredlist.org](https://www.iucnredlist.org/))                                                                                                                             |
| Fish silhouettes     | [PhyloPic](https://www.phylopic.org/), public domain / CC0 (credits in `output/phylopic_credits.csv`)                                                                                                                    |

## System requirements

R ≥ 4.1. Tested with R 4.3.2 on Windows 11. Required packages are listed in
`script/000_library.R` and installed automatically when missing. No
non-standard hardware is required.

## Reproducing the results

Open `FishUsage.Rproj` and run `script/00_RunAll.R` (or
`source("script/00_RunAll.R")` from the project root). The full pipeline runs in
about 15 minutes on a standard desktop computer. `script/00_make_figures.R`
redraws the figures only. Each numbered script can also be run on its own.

Steps that take hours (web scraping, missForest imputations, null models with
999 randomizations) are commented out and flagged `[LONG]`; their results are
stored in `dataPrepared/` and `output/` and reloaded instead.

Figures are written to `figures/` at their print size.

## Contact

Pierre Bouchet, CRBE, Université de Toulouse, France
(<pierre.bouchet@utoulouse.fr>, <pierrebdef@gmail.com>,
<https://pierrolaloune.github.io/>)
