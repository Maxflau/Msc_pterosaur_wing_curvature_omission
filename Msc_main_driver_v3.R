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
source("R/Msc_setup_01_libraries.R")
source("R/Msc_setup_02_biomechanics_functions.R")
# --- 2. Phylogenetic tree and specimen data ---
source("R/Msc_data_01_tree_creation.R")
source("R/Msc_data_02_specimen_table.R")
# --- 3. Wing outlines ---
source("R/Msc_outline_01_extract.R")
source("R/Msc_outline_02_load.R")
source("R/Msc_outline_03_orientation.R")
# --- 4. Flight performance metrics and theoretical spaces ---
source("R/Msc_performance_01_metrics.R")
source("R/Msc_performance_02_morphospace_pca.R")
source("R/Msc_performance_03_pareto_front.R")
source("R/Msc_performance_04_outline_efourier.R")
source("R/Msc_performance_05_theoretical_morphospace.R"); source("R/Msc_performance_06_theoretical_pareto_map.R")
# --- 5. General figures and visualizations ---
source("R/Msc_figures_01_setup.R")
source("R/Msc_figures_02_wing_background.R")
source("R/Msc_figures_03_colours_labels.R");source("R/Msc_figures_04_clade_colours_stress.R")
source("R/Msc_figures_05_biomechanical_panels.R"); source("R/Msc_figures_06_group_morphospace.R")
# --- 6. Phylogenetic analyses ---
source("R/Msc_phylo_01_prep.R"); source("R/Msc_phylo_02_signal_pgls.R"); source("R/Msc_phylo_03_morphospace_figures.R");
source("R/Msc_phylo_04_grade_comparison.R"); source("R/Msc_phylo_05_tree_figures.R")
# --- 7. Temporal (through-time) analyses ---
source("R/Msc_temporal_01_time_bins.R"); source("R/Msc_temporal_02_disparity_hulls.R"); source("R/Msc_temporal_03_summary_plots.R");source("R/Msc_temporal_04_multipanel_plots.R"); source("R/Msc_temporal_05_updated_plots.R"); source("R/Msc_temporal_06_size_disparity_panels.R")
# --- 8. Clade-level statistics and PCA ---
source("R/Msc_clade_01_save_tables.R");source("R/Msc_clade_02_summary.R");source("R/Msc_clade_03_pca_permanova.R")
source("R/Msc_stats_00_setup.R"); source("R/Msc_stats_01_group_metric_summary.R")
source("R/Msc_stats_02_kruskal_dunn.R"); source("R/Msc_stats_03_disparity_by_group.R");
# --- 9. Environment and depositional context (palaeoecology) ---
source("R/Msc_environment_01_summary.R"); source("R/Msc_environment_02_comparison.R"); source("R/Msc_environment_03_preservation_test.R")
# --- 10. Statistical tests (PERMANOVA, variance, allometry) ---
source("R/Msc_stats_04_permanova.R");source("R/Msc_stats_05_permanova_pc_scores.R");source("R/Msc_stats_06_permanova_pca_tables.R");
source("R/Msc_stats_07_pairwise_variance_loadings.R")
source("R/Msc_stats_08_allometry_data_prep.R"); source("R/Msc_stats_09_skeletal_allometry.R"); source("R/Msc_stats_10_wing_vs_skeleton.R")
source("R/Msc_stats_11_phylo_vs_nonphylo.R");source("R/Msc_stats_12_violin_plots.R")
source("R/Msc_stats_13_curvature_distribution.R")
source("R/Msc_stats_14_depositional_data_prep.R");source("R/Msc_stats_15_preservation_figures.R"); source("R/Msc_stats_16_taphonomic_sensitivity.R")
source("R/Msc_stats_17_optimality_through_time.R"); source("R/Msc_stats_18_optimality_bootstrap.R"); source("R/Msc_stats_19_categorical_associations.R")
# --- 11. Outline shape analysis (EFA) ---
source("R/Msc_efa_00_parse_nts.R")
source("R/Msc_efa_01_data_join.R"); source("R/Msc_efa_02_plot_helpers.R")
source("R/Msc_efa_03_base_theme.R"); source("R/Msc_efa_04_morphospace_figures.R")
source("R/Msc_efa_05_group_morphospace_figures.R"); source("R/Msc_efa_06_phylo_morphospace.R"); source("R/Msc_efa_07_grade_comparison_figures.R")
source("R/Msc_efa_08_permanova_shape.R"); source("R/Msc_efa_09_biomechanical_distance.R")
source("R/Msc_efa_10_distance_tests.R"); source("R/Msc_efa_11_specimen_table.R")
source("R/Msc_efa_12_save_helpers.R"); source("R/Msc_efa_13_loadings_variance_tables.R")
source("R/Msc_efa_14_permanova_all_groups.R"); source("R/Msc_efa_15_environment_curvature_tests.R")
source("R/Msc_efa_16_biomech_shape_regression.R"); source("R/Msc_efa_17_phylo_signal_violins.R")
source("R/Msc_efa_18_time_bin_regression.R"); source("R/Msc_efa_19_morphospace_expansion.R")
source("R/Msc_efa_20_pca_variance_loadings.R")
# --- 12. Final results (morphospace, wing classification, PGLS) ---
source("R/Msc_final_01_stacked_morphospace.R")
source("R/Msc_final_02_wing_classification.R")
source("R/Msc_final_03_pgls_tests.R")