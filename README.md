# Msc_pterosaur_wing_curvature_omission
This repository contains the full pipeline, metadata (both pterosaur wing outline and paleoecological informations) to run the analysis and reproduce the supplementary material 

## Structure

- `Msc_main_driver_v3.R` — main entry point of the analysis pipeline. Run this script to reproduce the full analysis; it sources every module in the order required by the pipeline.
- `R/` — all source modules (`source()`d by the main driver), covering wing outline extraction/loading, data handling, theoretical morphospaces, phylogenetics, and statistical analyses (PERMANOVA, allometry, disparity, EFA, PGLS, etc.).
- `data/` — input datasets used throughout the pipeline:
  - `pteros_main_data.csv` / `pteros_main_data_v2.csv` — main specimen/trait database.
  - `Pterosaur_Outlines_Clean.nts` — pre-extracted wing outline coordinates (used by `Msc_outline_loading.R`; regenerate from raw silhouette images with `Msc_outline_extract_1.R` if needed).
  - `pca_loadings.csv`, `clade_summary_statistics.csv`, `curvature_presence_summary.csv`, `disparity_by_environment.csv`, `allometry_correlation_output.csv`, `anatomical_wing_optimization_results.csv`, `performance_metrics_COMPLETE.csv`, `MarSoar.csv`, `Migrate.csv` — intermediate/supporting tables consumed and/or refreshed by the statistical modules.
  - `data/phylogenetics/` — phylogenetic inputs and outputs: `Henry` (base reference tree, NEXUS) and `Henry_updated.tre` / `Henry_updated.nex` (tree grafted with additional taxa, produced by `Msc_tree_creation.R` and consumed by the phylogenetic/PGLS/EFA modules).

Raw wing silhouette images used for outline extraction are not tracked in this repository (large binary assets); place them under `data/pteros_out_finale/` locally to re-run extraction from scratch, otherwise the pipeline uses the pre-extracted outlines in `data/Pterosaur_Outlines_Clean.nts`.

## Usage

1. Clone this repository.
2. Open `Msc_main_driver_v3.R` in RStudio (or set your working directory to the repository root).
3. Run the script — it will source all modules in `R/` in the correct order and execute the full analysis pipeline, reading its inputs from `data/`.
