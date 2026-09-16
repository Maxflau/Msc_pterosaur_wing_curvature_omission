# Msc_pterosaur_wing_curvature_omission
This repository contains the full pipeline, metadata (both pterosaur wing outline and paleoecological informations) to run the analysis and reproduce the supplementary material 

## Structure

- `Msc_main_driver_v3.R` — main entry point of the analysis pipeline. Run this script to reproduce the full analysis; it sources every module in the order required by the pipeline.
- `R/` — all source modules (`source()`d by the main driver), covering data handling, outline processing, theoretical morphospaces, phylogenetics, and statistical analyses (PERMANOVA, allometry, disparity, EFA, PGLS, etc.).

## Usage

1. Clone this repository.
2. Open `Msc_main_driver_v3.R` in RStudio (or set your working directory to the repository root).
3. Run the script — it will source all modules in `R/` in the correct order and execute the full analysis pipeline.
