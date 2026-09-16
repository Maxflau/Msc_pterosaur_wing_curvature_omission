WL_COL <- "Wing.loading..EWM."
results_list <- vector("list", length(outlines_list))
n_skipped <- 0
for (i in seq_along(outlines_list)) {
  
  specimen_name <- names(outlines_list)[i]
  outline <- outlines_list[[i]]
  
  specimen_row <- complete_data[which(complete_data$SPECIES..Pegas. == specimen_name), ]
  if (nrow(specimen_row) == 0) { n_skipped <- n_skipped + 1; next }
  specimen_row <- specimen_row[1, ]
  AR <- if ("WAR" %in% colnames(specimen_row)) as.numeric(specimen_row$WAR) else NA_real_
  WL <- if (WL_COL %in% colnames(specimen_row)) as.numeric(specimen_row[[WL_COL]]) else NA_real_
  MASS_KG   <- if ("Mass_kg" %in% colnames(specimen_row)) as.numeric(specimen_row$Mass_kg) else NA_real_
  WSPAN_CM  <- if ("Wingspan_cm" %in% colnames(specimen_row)) as.numeric(specimen_row$Wingspan_cm) else NA_real_
  r2_hat_val <- calculate_r2_hat(outline)
  metrics <- tryCatch({
    data.frame(
      species          = specimen_name,
      aspect_ratio     = AR,
      wing_loading     = WL,
      
      r2_hat           = r2_hat_val,
      
      # von_mises_stress is the ONLY stress metric in this pipeline.
      # stress_index/stress_index_DIAGNOSTIC and the old stress_root have
      # been removed entirely - no fallback, no diagnostic comparison.
      von_mises_stress = calculate_von_mises_stress(outline, MASS_KG, WSPAN_CM),
      
      wing_curvature   = calculate_wing_curvature(outline),
      shape_complexity = calculate_shape_complexity(outline),
      
      aspect_ratio_outline = calculate_aspect_ratio(outline),
      
      second_moment_DIAGNOSTIC = calculate_second_moment_DIAGNOSTIC(outline),
      
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(paste("Error calculating metrics for", specimen_name, ":", e$message, "\n"))
    NULL
  })
  
  results_list[[i]] <- metrics
  
  if (i %% 20 == 0) cat(paste("Processed", i, "/", length(outlines_list), "specimens\n"))
}

performance_data <- dplyr::bind_rows(results_list)

cat(paste("\nCalculated performance metrics for", nrow(performance_data), "specimens\n"))
cat(paste("Outlines with no matching database row:", n_skipped, "\n\n"))

if (nrow(performance_data) == 0 || !"species" %in% colnames(performance_data)) {
  stop("performance_data is empty or has no species column. Every metric call ",
       "failed - check the error lines printed above.")
}

# --------------------------------------------------------------------------------
# Validation - check these numbers before running anything downstream
# --------------------------------------------------------------------------------
cat("r2_hat (expected roughly 0.30-0.70; a spread over orders of magnitude\n")
cat("means outlines were not span-aligned):\n")
print(summary(performance_data$r2_hat))

cat("\nvon_mises_stress:\n")
print(summary(performance_data$von_mises_stress))

cat("\ntip_angle (degrees, expect roughly 50-120 per Walters et al.):\n")
print(summary(performance_data$tip_angle))

cat("\npitch_agility (r2_hat-based):\n")
print(summary(performance_data$pitch_agility))

wa <- complete_data$Wing_area_cm2[match(performance_data$species,
                                        complete_data$SPECIES..Pegas.)]
ok <- is.finite(performance_data$r2_hat) & is.finite(wa) & wa > 0
if (sum(ok) > 10) {
  cat(sprintf("\nCorrelation r2_hat vs real wing area (log-log): %.3f\n",
              cor(log(performance_data$r2_hat[ok]), log(wa[ok]))))
  cat("Near zero = scale-invariant, as intended.\n")
}

ok2 <- is.finite(performance_data$aspect_ratio) &
  is.finite(performance_data$aspect_ratio_outline)
if (sum(ok2) > 10) {
  cat(sprintf("\nCorrelation database AR vs outline AR: %.3f\n",
              cor(performance_data$aspect_ratio[ok2],
                  performance_data$aspect_ratio_outline[ok2])))
}
cat("\n")

# --------------------------------------------------------------------------------
# Merge with taxonomic and palaeoenvironmental data
# --------------------------------------------------------------------------------
join_cols <- c("SPECIES..Pegas.", "Order", "clade", "family", "Environment","Diet.1", "Diet.2", "Midpoint", "Period.Name", "Depositional.settings.paleoenvironment","Wingspan_cm", "Mass_kg", "Wing_area_cm2")
if ("Flight_category" %in% colnames(complete_data)) {
  join_cols <- c(join_cols, "Flight_category")
}
join_cols <- join_cols[join_cols %in% colnames(complete_data)]

performance_data <- performance_data %>%
left_join(complete_data %>% dplyr::select(all_of(join_cols)) %>% distinct(),by = c("species" = "SPECIES..Pegas."))

# --------------------------------------------------------------------------------
# Mass unit error detection and correction (g mistakenly entered as kg),
# using an external wingspan-based mass estimator rather than a regression
# fit on this (partly flawed) dataset itself:
#   Mass_est_kg = 0.2889 * Wingspan_m^2.6562
# Verified against two known specimens: Quetzalcoatlus (9.7 m) -> ~121 kg,
# matching published estimates (~120 kg); Dimorphodon (1.337 m) -> ~0.62 kg,
# matching published estimates (~1-2 kg) and confirming the recorded 625 in
# this dataset was a g-as-kg error. A specimen is corrected only if dividing
# its recorded Mass_kg by 1000 brings it close to this external estimate -
# nothing is divided blindly.
# --------------------------------------------------------------------------------
fix_mass_unit_errors <- function(data, mass_col = "Mass_kg", span_col = "Wingspan_cm",
                                 ratio_flag = 20, ratio_fix_tol = 5) {
  
  ok <- is.finite(data[[mass_col]]) & data[[mass_col]] > 0 &
    is.finite(data[[span_col]]) & data[[span_col]] > 0
  
  wingspan_m <- data[[span_col]] / 100
  mass_est   <- 0.2889 * wingspan_m^2.6562
  
  ratio <- data[[mass_col]] / mass_est
  flagged <- ok & is.finite(ratio) & ratio > ratio_flag
  
  corrected_mass <- data[[mass_col]]
  fixed <- rep(FALSE, nrow(data))
  
  for (i in which(flagged)) {
    candidate <- data[[mass_col]][i] / 1000
    candidate_ratio <- candidate / mass_est[i]
    if (candidate_ratio >= 1 / ratio_fix_tol && candidate_ratio <= ratio_fix_tol) {
      corrected_mass[i] <- candidate
      fixed[i] <- TRUE
    }
  }
  
  data[[paste0(mass_col, "_corrected")]] <- corrected_mass
  data$mass_est_kg         <- mass_est
  data$mass_unit_flagged   <- flagged
  data$mass_unit_corrected <- fixed
  
  cat(sprintf("\nMass unit check (external wingspan-based estimator): %d specimens flagged (ratio > %gx)\n",
              sum(flagged), ratio_flag))
  cat(sprintf("  -> %d corrected (/1000 brought them within %gx of the estimate)\n",
              sum(fixed), ratio_fix_tol))
  cat(sprintf("  -> %d still flagged, NOT corrected - inspect manually\n",
              sum(flagged) - sum(fixed)))
  
  if (sum(flagged) > 0) {
    cat("\nFlagged specimens:\n")
    print(data.frame(
      species          = data$species[flagged],
      Mass_kg_original  = round(data[[mass_col]][flagged], 3),
      Mass_kg_estimated = round(mass_est[flagged], 3),
      Mass_kg_corrected = round(corrected_mass[flagged], 3),
      Wingspan_cm       = data[[span_col]][flagged],
      corrected         = fixed[flagged]
    ))
  }
  cat("\n")
  
  data
}

if (all(c("Mass_kg", "Wingspan_cm") %in% colnames(performance_data))) {
  performance_data <- fix_mass_unit_errors(performance_data)
} else {
  cat("Mass_kg/Wingspan_cm absent - mass unit check skipped.\n")
  performance_data$Mass_kg_corrected <- performance_data$Mass_kg
}
# --------------------------------------------------------------------------------
# Wing loading residual and ratio - allometric size correction
# (log WL ~ log Mass_kg). wing_loading_ratio = exp(residual): positive,
# usable directly in the log-transformed PCA further down the pipeline.
# --------------------------------------------------------------------------------
if (all(c("wing_loading", "Mass_kg") %in% colnames(performance_data))) {
  ok_wl <- is.finite(performance_data$wing_loading) & performance_data$wing_loading > 0 &
    is.finite(performance_data$Mass_kg) & performance_data$Mass_kg > 0
  if (sum(ok_wl) > 10) {
    m_wl <- lm(log(wing_loading) ~ log(Mass_kg), data = performance_data[ok_wl, ])
    
    performance_data$wing_loading_resid <- NA_real_
    performance_data$wing_loading_resid[ok_wl] <- residuals(m_wl)
    
    performance_data$wing_loading_ratio <- NA_real_
    performance_data$wing_loading_ratio[ok_wl] <- exp(residuals(m_wl))
    
    cat("\nwing_loading_resid (additive, descriptive only, NOT in the PCA):\n")
    print(summary(performance_data$wing_loading_resid))
    
    cat("\nwing_loading_ratio (multiplicative, positive, usable in the PCA):\n")
    print(summary(performance_data$wing_loading_ratio))
    cat(sprintf("Correlation with log(Mass_kg): %.4f (should be ~0)\n",
                cor(log(performance_data$wing_loading_ratio[ok_wl]),
                    log(performance_data$Mass_kg[ok_wl]))))
  } else {
    cat("\nToo few finite wing_loading/Mass_kg pairs - wing_loading_ratio not computed.\n")
  }
} else {
  cat("\nwing_loading or Mass_kg absent - wing_loading_ratio not computed.\n")
}

cat(paste("\nFinal dataset:", nrow(performance_data), "rows x",
          ncol(performance_data), "columns\n\n"))