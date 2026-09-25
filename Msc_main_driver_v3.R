# ============================================================
# MAIN SCRIPT - runs the whole analysis in one go.
# Just run this file: each "source(...)" line below loads and
# runs one pipeline module, in the exact order required
# (do not change this order).
# ============================================================
rm(list=ls()) # clear R's memory before starting
# Automatically move to the repository folder (if opened in RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
# --- 1. Basic tools: packages and shared functions ---
source("R/Msc_librairies.R")
source("R/Msc_biomechanical_function_v3.R")
# --- 2. Phylogenetic tree and specimen data ---
source("R/Msc_tree_creation.R")
source("R/Msc_data_handling.R")
# --- 3. Wing outlines ---
source("R/Msc_outline_extract_1.R")
source("R/Msc_outline_loading.R")
source("R/Msc_orien_v1.R")
# --- 4. Flight performance metrics and theoretical spaces ---
source("R/Msc_performance_metrics_v3.R")
source("R/Msc_impossible_regions_v3.R")
source("R/Msc_pareto_front_1.5.R")
source("R/Msc_rep.R")
source("R/Msc_theoretical_spaces_v3.R"); source("R/Msc_theoretical_space_v3.2.R")
# --- 5. General figures and visualizations ---
source("R/Msc_plot_layer_1.R")
source("R/Msc_data_visual_v3.R")
source("R/Msc_fig_prep_1.R");source("R/Msc_fig_prep_2.R")
source("R/Msc_biomech_plot.R"); source("R/Msc_aBio_plot.R")
# --- 6. Phylogenetic analyses ---
source("R/Msc_phylo_v3.1.R"); source("R/Msc_phylo_v3.2.R"); source("R/Msc_phylo_v3.3.R");
source("R/Msc_phylo_v3.4.R"); source("R/Msc_phylo_front_1.R"); source("R/Msc_phylo_signal_v3.R")
# --- 7. Temporal (through-time) analyses ---
source("R/Msc_temp_3.1.R"); source("R/Msc_temp_3.1.5.R"); source("R/Msc_temp_3.2.R");source("R/Msc_temp_3.3.R"); source("R/Msc_temp_3.4.R"); source("R/Msc_temp_3.4_5.R")
# --- 8. Clade-level statistics and PCA ---
source("R/Msc_stats_2.R");source("R/Msc_clade_sum_2.R");source("R/Msc_PCA_loading_2.R")
source("R/Msc_stat_prep.R"); source("R/Msc_stat_metric.R")
source("R/Msc_stat_KW_Dunn.R"); source("R/Msc_stat_disparity.R");
# --- 9. Environment and depositional context (palaeoecology) ---
source("R/Msc_code_environment_control.R"); source("R/Msc_code_environment_control_2.R"); source("R/Msc_Depo_stats.R")
# --- 10. Statistical tests (PERMANOVA, variance, allometry) ---
source("R/Msc_stat_permanova.R");source("R/Msc_stat_permanova_2.R");source("R/Msc_permanova_3.R");
source("R/Msc_stat_variance_loading.R")
source("R/Msc_stat_allometry_1.R"); source("R/Msc_stat_allometry_2.R"); source("R/Msc_stat_allometry_3.R")
source("R/Msc_stat_checking.R");source("R/Msc_stat_violin.R")
source("R/Msc_curvature_distribution_v2.R")
source("R/Msc_stat_Depo_1.R");source("R/Msc_stat_Depo_2.R"); source("R/Msc_stat_Depo_3.R")
source("R/Msc_optimality_temp_1.R"); source("R/Msc_optimality_temp_2.R"); source("R/Msc_heatmap_stat.R")
# --- 11. Outline shape analysis (EFA) ---
source("R/Msc_parse_file.R")
source("R/MSc_EFA_2.1.R"); source("R/MSc_EFA_2.2.R")
source("R/MSc_EFA_2.3.R"); source("R/MSc_EFA_2.4.R")
source("R/MSc_EFA_2.5.R"); source("R/MSc_EFA_2.6.R"); source("R/MSc_EFA_2.6.5.R")
source("R/MSc_EFA_2.7.R"); source("R/MSc_EFA_2.8.R")
source("R/MSc_EFA_2.9.R"); source("R/MSc_EFA_2.10.R")
source("R/MSc_EFA_2.11.R"); source("R/MSc_EFA_2.12.R")
source("R/MSc_EFA_2.13.R"); source("R/MSc_EFA_2.14.R")
source("R/MSc_EFA_2.15.R"); source("R/MSc_EFA_2.16.R")
source("R/MSc_EFA_2.17.R"); source("R/MSc_EFA_2.18.R")
# --- 12. Final results (morphospace, wing classification, PGLS) ---
source("R/Msc_stacked_morphospace.R")
source("R/Msc_wing_classification.R")
source("R/Msc_PGLS_2.R")