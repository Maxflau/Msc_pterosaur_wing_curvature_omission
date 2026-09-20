# Build time-bin summaries and tracked variables for the temporal figures (Part A).
cat("\n================================================================================\n")
cat("TEMPORAL METRICS - PART A (time bins, tracked variables, bootstrap)\n")
cat("================================================================================\n\n")

N_BOOT <- 2000
set.seed(2026)

#####################################################################################
# 1. Put each specimen into the broad time bins used in the figures
#####################################################################################
temporal_data <- performance_data_clean %>%
  dplyr::filter(!is.na(Midpoint), !is.na(Period.Name), !is.na(clade))

temporal_data$Period_Combined <- temporal_data$Period.Name
temporal_data$Period_Combined[temporal_data$Period_Combined %in%
                                c("Early Jurassic", "Middle Jurassic")] <-"Early-Middle Jurassic"
temporal_data$Time_Bin <- factor( temporal_data$Period_Combined,
                                  levels = c("Late Triassic", "Early-Middle Jurassic", "Late Jurassic",
                                             "Early Cretaceous", "Late Cretaceous"))
temporal_data <- temporal_data %>% dplyr::filter(!is.na(Time_Bin))

cat(sprintf("Specimens with temporal data: %d\n", nrow(temporal_data)))
print(table(temporal_data$Time_Bin))
cat("\nNOTE: unbalanced bins. Interpret bin means with their intervals, never alone.\n\n")

######################################################################################
# 2. List every variable that will be followed through time
######################################################################################
SHAPE_VARS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress",
                          "wing_curvature", "shape_complexity",
                          "wing_loading_ratio"), colnames(temporal_data))

SCORE_VARS <- intersect(c("pareto_score_aspect_ratio", "pareto_score_r2_hat",
                          "pareto_score_von_mises_stress",
                          "pareto_score_wing_loading_ratio"), colnames(temporal_data))

TRADEOFF_VARS <- intersect(c("pareto_tradeoff_aspect_ratio", "pareto_tradeoff_r2_hat",
                             "pareto_tradeoff_von_mises_stress",
                             "pareto_tradeoff_wing_loading_ratio"), colnames(temporal_data))

SIZE_VARS  <- intersect(c("Wingspan_cm", "Mass_kg", "Wing_area_cm2"),
                        colnames(temporal_data))

OPT_VARS   <- intersect("pareto_rank_ratio", colnames(temporal_data))

TRACKED <- c(SHAPE_VARS, SCORE_VARS, TRADEOFF_VARS, SIZE_VARS, OPT_VARS)

if (length(SIZE_VARS) == 0) {
  warning("No size variables found. Add Wingspan_cm, Mass_kg and Wing_area_cm2 ",
          "to join_cols in Msc_performance_metrics_v3.R.")
}
if (length(SCORE_VARS) == 0) {
  cat("NOTE: pareto_score_* columns absent - source Msc_pareto_front_1.R ",
      "before this script.\n", sep = "")
}
cat("Tracked variables:\n")
cat("  shape:", paste(SHAPE_VARS, collapse = ", "), "\n")
cat("  pareto score:", paste(SCORE_VARS, collapse = ", "), "\n")
cat("  pareto trade-off:", paste(TRADEOFF_VARS, collapse = ", "), "\n")
cat("  size: ", paste(SIZE_VARS, collapse = ", "), "\n")
cat("  optimality:", paste(OPT_VARS, collapse = ", "), "\n\n")

#####################################################################################
# 3. Calculate a mean and uncertainty range for each time bin
#####################################################################################
boot_mean_ci <- function(x, n_boot = N_BOOT) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(mean = NA, lo = NA, hi = NA, n = length(x)))
  bm <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), lo = unname(quantile(bm, 0.025)),
    hi = unname(quantile(bm, 0.975)), n = length(x))
}

temporal_metrics <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
  d <- temporal_data[temporal_data$Time_Bin == tb, ]
  if (nrow(d) == 0) return(NULL)

  # Summarise size on a log scale so very large values do not dominate the mean.
  do.call(rbind, lapply(TRACKED, function(v) {
    x <- d[[v]]
    if (v %in% SIZE_VARS) x <- log10(x[is.finite(x) & x > 0])
    ci <- boot_mean_ci(x)
    data.frame(Time_Bin = tb, variable = v, scale = if (v %in% SIZE_VARS) "log10" else "raw",
               mean = ci["mean"], lo = ci["lo"], hi = ci["hi"], n = ci["n"],
               stringsAsFactors = FALSE)
  }))
}))

temporal_metrics$Time_Bin <- factor(temporal_metrics$Time_Bin,
                                    levels = levels(temporal_data$Time_Bin))
rownames(temporal_metrics) <- NULL

cat("Bootstrapped means by time bin:\n")
print(temporal_metrics, digits = 4)
cat("\nAvailable: temporal_data, temporal_metrics, N_BOOT\n\n")
