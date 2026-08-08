# FishUsage

[![Code DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21848876.svg)](https://doi.org/10.5281/zenodo.21848876)

Data, scripts and workflows for a study of human uses of freshwater fishes.

The analysis covers **8,970 freshwater fish species** and five categories of human
use (fisheries, aquaculture, aquarium trade, recreational fishing). It asks
where used species sit in the morphological space of freshwater fishes, whether
morphologically distinctive species are more likely to be targeted, and how much
of that space would be lost if threatened species disappeared.

---

## Repository structure

```text
FishUsage/
├── script/                   # Analysis pipeline, run in numeric order
│   ├── 000_library.R                 # Packages: installs what is missing, then loads
│   ├── 000_functions.R               # Every custom function, no analysis of its own
│   ├── 000_LoadDataR.R               # Traits, phylogeny, IUCN status, human uses
│   ├── 000_ScrappingData.R           # Human uses scraped from FishBase pages
│   ├── 01_*.R … 07_*.R               # Analyses: functional richness, null models,
│   │                                 #   distinctiveness, imputation error
│   ├── 08_*.R … 12_*.R               # Figures: functional spaces, richness loss,
│   │                                 #   distinctiveness, PCA loadings, deficit maps
│   ├── 13_*.R, 100_*.R               # Sensitivity to the trait imputation
│   ├── plot_pca_correlation_circle.R # Alternative correlation circles
│   ├── web scrapping percentage.R    # What the scraping added over rfishbase
│   └── test_new_funspace.R           # Exploratory, not part of the final analysis
│
├── dataOriginal/             # Source datasets, as downloaded
│                             #   FISHMORPH traits and phylogeny, IUCN Red List
│
├── dataPrepared/Fish/        # Intermediate tables built by 000_LoadDataR.R
│                             #   cleaned traits, imputed traits, phylogenetic PCoA
│
├── output/                   # Precomputed results reloaded by the scripts (~150 MB)
│
├── figures/                  # Main and supplementary figures
│   ├── Clean/                # Final versions
│   └── individual_panels/    # One square panel per use
│
├── README.md                 # This file
└── .gitignore
```

Every script uses paths relative to the project root, so R must be started from
the folder that contains `script/`, `dataOriginal/`, `dataPrepared/`, `figures/`
and `output/`.

---

## Data availability

The datasets and all precomputed results are included in this repository, under
`dataOriginal/`, `dataPrepared/` and `output/`. They are also archived on Zenodo
together with the code: [10.5281/zenodo.21848876](https://doi.org/10.5281/zenodo.21848876).

Original sources:

| Dataset | Source |
| --- | --- |
| Morphological traits and phylogeny | Brosse, S. et al. FISHMORPH: A global database on morphological traits of freshwater fishes. *Global Ecol. Biogeogr.* **30**, 2330–2336 (2021). |
| Length, weight, human uses | Froese, R. & Pauly, D. FishBase. (2025). Accessed through [`rfishbase`](https://docs.ropensci.org/rfishbase/) and by scraping the summary pages. |
| Conservation status | IUCN. The IUCN Red List of Threatened Species. (2025). |

The datasets redistributed under `dataOriginal/` remain subject to the terms of
their respective providers. Please consult those terms before reusing or
redistributing them.

---

## Reproducing the analysis

Requires **R ≥ 4.1** (the native pipe `|>` is used in places).

1. Clone the repository:

   ```bash
   git clone https://github.com/pierrolaloune/FishUsage.git
   cd FishUsage
   ```

2. Open R from the project root, then source the setup scripts in this order, at
   the start of every session:

   ```r
   source("script/000_library.R")    # installs any missing package, then loads all of them
   source("script/000_functions.R")  # defines every custom function
   source("script/000_LoadDataR.R")  # builds the trait, IUCN and human-use tables
   ```

3. Run any numbered script. Each one reloads what it needs from `output/`, so
   they are independent and can be run in any order:

   ```r
   source("script/01_FRic_Dissim.R")
   ```

   The one exception is `web scrapping percentage.R`, which reuses an object
   built by `000_LoadDataR.R` and must run in the same session.

Figures are written to `figures/`, result tables to `output/`.

> **A note on runtime.** Web scraping, random-forest imputation and the null
> models with 999 replicates each take hours. Every one of those steps is
> **already commented out**, with its result stored in `output/` or
> `dataPrepared/` and reloaded on the next line, so the scripts run end to end as
> they are. Each is flagged `[LONG]` in the section title of the script it
> belongs to. Uncomment a block only to rebuild that file from scratch.

---

## How to cite

If you use this code or these data, please cite the archive:

> Bouchet, P. pierrolaloune/FishUsage. Zenodo. https://doi.org/10.5281/zenodo.21848876

---

## Contact

Pierre Bouchet, CRBE, Université de Toulouse, France

- Email: pierre.bouchet@utoulouse.fr
- Website: <https://pierrolaloune.github.io/>
