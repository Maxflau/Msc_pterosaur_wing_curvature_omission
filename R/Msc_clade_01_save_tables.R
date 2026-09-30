# Save the main result tables that later scripts may reuse.
# Save the full cleaned table in the supplement folder.
write.csv(performance_data_clean, "output/results/Supplemental_performance/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("✓ Exported: output/results/performance_metrics_complete_CLADE.csv\n")

# Save the same full cleaned table in the main results folder.
write.csv(performance_data_clean, "output/results/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("✓ Exported: output/results/performance_metrics_complete_CLADE.csv\n")

# Save the disparity table when an earlier script has already created it.
if (exists("disparity_metrics")) {
  write.csv(disparity_metrics, 
            "output/results/Supplemental_performance/disparity_metrics.csv", row.names = FALSE)
  cat("✓ Exported: output/results/disparity_metrics.csv\n")
} else {
  cat("SKIPPED disparity_metrics - not found (source Msc_temp_3_1.R first).\n")
}

# Save the time-bin summary table when it is available.
if (exists("temporal_metrics")) {
  write.csv(temporal_metrics, 
            "output/results/Supplemental_performance/temporal_metrics_CLADE.csv", row.names = FALSE)
  cat("✓ Exported: output/results/Supplemental_performance/temporal_metrics_CLADE.csv\n")
} else {
  cat("SKIPPED temporal_metrics - not found (source Msc_temp_3_1.R first).\n")
}
# Save the theoretical shape table when it is available.
if (exists("theoretical_data")) {
  write.csv(theoretical_data, 
            "output/results/Supplemental_performance/theoretical_shape_space.csv", row.names = FALSE)
  cat("✓ Exported: output/results/theoretical_shape_space.csv\n\n")
} else {
  cat("SKIPPED theoretical_data - not found (source Msc_theoretical_spaces_v3.R first).\n\n")
}
