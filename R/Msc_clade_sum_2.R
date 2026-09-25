# Summarise the main biomechanical measurements for each clade.
# These tables are saved for later reporting.

# --- 1. Mean and spread by clade ----------------------------------------------
summary_by_clade <- performance_data_clean %>%
  group_by(clade) %>%
  summarise(
    n = n(),
    mean_AR = mean(aspect_ratio, na.rm = TRUE),
    sd_AR = sd(aspect_ratio, na.rm = TRUE),
    mean_r2hat = mean(r2_hat, na.rm = TRUE),
    sd_r2hat = sd(r2_hat, na.rm = TRUE),
    mean_WL = mean(wing_loading, na.rm = TRUE),
    sd_WL = sd(wing_loading, na.rm = TRUE),
    mean_WL_ratio = mean(wing_loading_ratio, na.rm = TRUE),
    sd_WL_ratio = sd(wing_loading_ratio, na.rm = TRUE),
    mean_stress = mean(von_mises_stress, na.rm = TRUE),
    sd_stress = sd(von_mises_stress, na.rm = TRUE),
    mean_curvature = mean(wing_curvature, na.rm = TRUE),
    sd_curvature = sd(wing_curvature, na.rm = TRUE),
    mean_complexity = mean(shape_complexity, na.rm = TRUE),
    sd_complexity = sd(shape_complexity, na.rm = TRUE)
    # Add more measurements here only if later scripts also use them.
  )
write.csv(summary_by_clade,
          "output/results/Supplemental_performance/summary_statistics_by_clade.csv", row.names = FALSE)
cat("✓ Exported: output/results/summary_statistics_by_clade.csv\n\n")
# --- 2. Mean and spread by order ----------------------------------------------
summary_by_order <- performance_data_clean %>%
  group_by(Order) %>%
  summarise(
    n = n(),
    mean_AR = mean(aspect_ratio, na.rm = TRUE),
    sd_AR = sd(aspect_ratio, na.rm = TRUE),
    mean_r2hat = mean(r2_hat, na.rm = TRUE),
    sd_r2hat = sd(r2_hat, na.rm = TRUE),
    mean_WL = mean(wing_loading, na.rm = TRUE),
    sd_WL = sd(wing_loading, na.rm = TRUE),
    mean_WL_ratio = mean(wing_loading_ratio, na.rm = TRUE),
    sd_WL_ratio = sd(wing_loading_ratio, na.rm = TRUE),
    mean_stress = mean(von_mises_stress, na.rm = TRUE),
    sd_stress = sd(von_mises_stress, na.rm = TRUE),
    mean_curvature = mean(wing_curvature, na.rm = TRUE),
    sd_curvature = sd(wing_curvature, na.rm = TRUE),
    mean_complexity = mean(shape_complexity, na.rm = TRUE),
    sd_complexity = sd(shape_complexity, na.rm = TRUE)
  )
write.csv(summary_by_order, "output/results/Supplemental_performance/summary_statistics_by_order.csv", row.names = FALSE)
cat("✓ Exported: output/results/summary_statistics_by_order.csv\n\n")
