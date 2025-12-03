# FishUsage
### Data and code to reproduce the analyses from  
**“Targeting morphologically unique fishes amplifies the risk of functional erosion”**  
*Pierre Bouchet, Sébastien Brosse & Aurèle Toussaint*  
CNRS, Université de Toulouse

---

## 📖 Overview

This repository contains all data and R scripts necessary to reproduce the analyses presented in the manuscript *Targeting morphologically unique fishes amplifies the risk of functional erosion* (paper not yet published).

The project evaluates whether human uses of freshwater fishes — including food, ornamental trade, and recreational activities — are concentrated on narrow portions of the global fish functional spectrum. We test the expectation that humans select a non-random subset of species, and that morphologically unique species, often already threatened, are disproportionately targeted, increasing risks of functional erosion.

The repository covers **all freshwater fish species globally**, integrating phylogeny, functional morphology, human use information, and conservation status (IUCN 2024 assessments).

---

## 🗂️ Repository Structure

### **Data**
The repository includes both raw and processed datasets:

- **Phylogeny** (FishMorph)
- **Morphological traits** (FishMorph + FishBase)
- **IUCN Red List 2024 assessments**
- **Human use categories** (food, ornament, recreation)
- Formats: **CSV** and **RDS**

All datasets used here are open-access but **require proper citation**:
- **IUCN (2024). The IUCN Red List of Threatened Species. Version 2024-1. https://www.iucnredlist.org**
- **Froese, R. & Pauly, D. (eds.) (2024). FishBase. World Wide Web electronic publication. www.fishbase.org Accessible at: https://www.fishbase.org**
- **Brosse, S., Charpin N., Su G., Toussaint A., Herrera-R G. A., Tedesco P. A., & Villéger S. (2021). FISHMORPH: A global database on morphological traits of freshwater fishes. Global Ecology and Biogeography, 30, 2330–2336. https://doi.org/10.1111/geb.13395**

---

## 🧠 Scripts

Scripts follow a **numbered sequential pipeline** to allow full reproducibility:

- **000_library.R** — loads all required R packages  
- **000_functions.R** — functions used across the global analysis  
- **000_LoadData.R** — imports, cleans, and prepares all datasets  
- **000_ScrappingData.R** — optional script to reproduce the scraping workflow (not necessary to rerun)  
- **0XX_analysis_*.R** — numbered scripts performing the full analysis pipeline  
- **1XX_figures_*.R** — numbered scripts generating the manuscript and supplementary figures  

Some analyses are computationally heavy.  
**Long HPC-dependent steps are commented out**, and corresponding processed `.rds` files are provided so users can run the workflow without recomputing large objects.

---

## 🔄 Reproducibility Workflow (simple)

1. Clone the repository  
2. Open R (≥ **4.3.2**)  
3. Run `000_library.R`  
4. Run `000_functions.R`  
5. Run `000_LoadData.R`  
6. Run the numbered analysis scripts in order  
7. Run the figure scripts to reproduce the paper’s figures and supplementary materials  

---

## 📊 Figures

Scripts will reproduce:
- All **main figures** included in the paper  
- All **supplementary figures**  

No additional outputs are included in this repository.

---

## 📦 Requirements

### **R version**
This project was developed using:

