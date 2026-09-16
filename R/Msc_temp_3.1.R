cat("\n================================================================================\n")
cat("TEMPORAL METRICS\n")
cat("================================================================================\n\n")

N_BOOT <- 2000
set.seed(2026)

#####################################################################################
# 1. Time bins
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
# 2. Variables tracked through time
######################################################################################
SHAPE_VARS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress",
                          "wing_curvature", "shape_complexity",
                          "wing_loading_ratio"), colnames(temporal_data))

# Extended to all 4 Pareto objectives, matching Msc_pareto_front_1.R.
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
# 3. Bootstrap mean and percentile interval
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
  
  do.call(rbind, lapply(TRACKED, function(v) {
    # Size variables span orders of magnitude: summarise on the log scale
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
cat("\n")

######################################################################################
# 4. Disparity through time
######################################################################################
# Sum of variance and mean pairwise distance on the performance PC plane,
# bootstrapped as in Liu et al. Both are size-free by construction.
disparity_boot <- function(d, n_boot = 1000) {
  M <- cbind(d$PC1, d$PC2)
  M <- M[complete.cases(M), , drop = FALSE]
  if (nrow(M) < 4) return(c(sov = NA, sov_lo = NA, sov_hi = NA,
                            mpd = NA, mpd_lo = NA, mpd_hi = NA))
  sov <- function(m) sum(apply(m, 2, var))
  mpd <- function(m) mean(dist(m))
  
  bs <- replicate(n_boot, {
    idx <- sample(nrow(M), nrow(M), replace = TRUE)
    c(sov(M[idx, , drop = FALSE]), mpd(M[idx, , drop = FALSE]))
  })
  c(sov = sov(M), sov_lo = quantile(bs[1, ], 0.025), sov_hi = quantile(bs[1, ], 0.975),
    mpd = mpd(M), mpd_lo = quantile(bs[2, ], 0.025), mpd_hi = quantile(bs[2, ], 0.975))
}

disparity_metrics <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
  d <- temporal_data[temporal_data$Time_Bin == tb, ]
  if (nrow(d) == 0) return(NULL)
  cbind(data.frame(Time_Bin = tb, n = nrow(d)),
        as.data.frame(t(disparity_boot(d))))
}))
disparity_metrics$Time_Bin <- factor(disparity_metrics$Time_Bin,
                                     levels = levels(temporal_data$Time_Bin))
rownames(disparity_metrics) <- NULL

cat("Disparity by time bin (performance PC1-PC2 plane):\n")
print(disparity_metrics, digits = 4)
cat("\n")

######################################################################################
# 5. Convex hulls per clade and bin, for the occupation figure
######################################################################################
temporal_hulls_clade <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
  d <- temporal_data[temporal_data$Time_Bin == tb, ]
  do.call(rbind, lapply(unique(d$clade), function(cl) {
    dc <- d[d$clade == cl & is.finite(d$PC1) & is.finite(d$PC2), ]
    if (nrow(dc) < 3) return(NULL)
    h <- chull(dc$PC1, dc$PC2)
    data.frame(Time_Bin = tb, clade = cl, PC1 = dc$PC1[h], PC2 = dc$PC2[h])
  }))
}))
temporal_hulls_clade$Time_Bin <- factor(temporal_hulls_clade$Time_Bin, levels = levels(temporal_data$Time_Bin))
cat(sprintf("Clade hulls built: %d clade-bin combinations\n\n",nrow(unique(temporal_hulls_clade[, c("Time_Bin", "clade")]))))

cat("Available: temporal_data, temporal_metrics, disparity_metrics, temporal_hulls_clade\n\n")