write.csv(performance_data_clean, "output/results/supplementals/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("✓ Exported: output/results/performance_metrics_complete_CLADE.csv\n")

write.csv(performance_data_clean, "output/results/performance_metrics_complete_CLADE.csv", row.names = FALSE)
cat("✓ Exported: output/results/performance_metrics_complete_CLADE.csv\n")

# disparity_by_clade / temporal_by_clade / clade_time_verification do not
if (exists("disparity_metrics")) {
  write.csv(disparity_metrics, 
            "output/results/supplementals/disparity_metrics.csv", row.names = FALSE)
  cat("✓ Exported: output/results/disparity_metrics.csv\n")
} else {
  cat("SKIPPED disparity_metrics - not found (source Msc_temp_3_1.R first).\n")
}

# temporal_metrics_comprehensive -> temporal_metrics
if (exists("temporal_metrics")) {
  write.csv(temporal_metrics, 
            "output/results/supplementals/temporal_metrics_CLADE.csv", row.names = FALSE)
  cat("✓ Exported: output/results/supplementals/temporal_metrics_CLADE.csv\n")
} else {
  cat("SKIPPED temporal_metrics - not found (source Msc_temp_3_1.R first).\n")
}
# theoretical_shapes -> theoretical_data (built in Msc_theoretical_spaces_v3.R)
if (exists("theoretical_data")) {
  write.csv(theoretical_data, 
            "output/results/supplementals/theoretical_shape_space.csv", row.names = FALSE)
  cat("✓ Exported: output/results/theoretical_shape_space.csv\n\n")
} else {
  cat("SKIPPED theoretical_data - not found (source Msc_theoretical_spaces_v3.R first).\n\n")
}