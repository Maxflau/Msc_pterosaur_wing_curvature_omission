cat("\n--- S14/S15: WING METRICS AGAINST THE SKELETON ---\n\n")

# Test whether the wing measures still change with size or limb proportions.
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_06a_postcranial.R first.")

MIN_FIT <- 10

TARGETS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress", "wing_curvature",
                       "shape_complexity", "froude"),
                     colnames(POSTCRANIAL))

# Dimensionless ratios should not scale with size: expected slope 0.
# Wing loading is force per area and is expected to rise with size.
EXPECTED <- c(aspect_ratio = 0, r2_hat = 0, von_mises_stress = 0,
              wing_curvature = 0, shape_complexity = 0, froude = NA)

cat("Targets:", paste(TARGETS, collapse = ", "), "\n\n")

# --------------------------------------------------------------------------------
# 1. Metric against skeletal size
# --------------------------------------------------------------------------------
# Test each wing measure against several body-size proxies.
fit_metric <- function(target, size_var, d) {
  x <- d[[size_var]]; y <- d[[target]]
  ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  if (sum(ok) < MIN_FIT) return(NULL)

  m  <- lm(log(y[ok]) ~ log(x[ok]))
  cf <- summary(m)$coefficients
  slope <- cf[2, 1]; se <- cf[2, 2]
  exp_sl <- EXPECTED[[target]]

  p_iso <- if (is.na(exp_sl)) NA else
    2 * pt(abs((slope - exp_sl) / se), df = sum(ok) - 2, lower.tail = FALSE)

  data.frame(target = target, predictor = size_var, n = sum(ok),
             slope = round(slope, 4), se = round(se, 4),
             r_squared = round(summary(m)$r.squared, 4),
             expected_slope = exp_sl,
             p_vs_expected = if (is.na(p_iso)) NA else signif(p_iso, 4),
             stringsAsFactors = FALSE)
}

SIZE_PREDICTORS <- intersect(c("SkeletalSize", "WingLength", "HindLimb",
                               "Wingspan_cm", "Mass_kg"),
                             colnames(POSTCRANIAL))
cat("Size predictors:", paste(SIZE_PREDICTORS, collapse = ", "), "\n")

metric_allo <- do.call(rbind, lapply(TARGETS, function(tg)
  do.call(rbind, lapply(SIZE_PREDICTORS, function(sv)
    fit_metric(tg, sv, POSTCRANIAL)))))

if (is.null(metric_allo)) stop("No metric fit succeeded.")

metric_allo$p_adjusted_BH <- signif(p.adjust(metric_allo$p_vs_expected, "BH"), 4)
metric_allo$size_dependent <- !is.na(metric_allo$p_adjusted_BH) &
  metric_allo$p_adjusted_BH < 0.05

cat("\nWing metrics against size:\n")
# Show a fixed number of decimals so every row lines up.
metric_allo_display <- metric_allo
metric_allo_display$slope         <- format_fixed(metric_allo$slope, 4)
metric_allo_display$r_squared     <- format_fixed(metric_allo$r_squared, 4)
metric_allo_display$p_adjusted_BH <- format_fixed(metric_allo$p_adjusted_BH, 4)
print(metric_allo_display[, c("target", "predictor", "n", "slope", "r_squared",
                      "expected_slope", "p_adjusted_BH", "size_dependent")],
      row.names = FALSE)
write_supp(metric_allo_display, "S14_metric_allometry")

n_dep <- sum(metric_allo$size_dependent[metric_allo$expected_slope == 0], na.rm = TRUE)
cat(sprintf("\n%d of the dimensionless metrics still scale with size.\n", n_dep))
cat("Zero would confirm the non-dimensionalisation; a large number means a\n")
cat("metric retains a size component and should be reported as such.\n\n")

# --------------------------------------------------------------------------------
# 2. Metric against skeletal proportions
# --------------------------------------------------------------------------------
# Turn each bone into a size-free proportion, then test those against the metrics.
# Each element divided by the geometric mean: a size-free shape variable. A
# relationship here is a proportion effect, independent of how big the animal is.
PROP_ELEMENTS <- intersect(c("Humerus", "Ulna", "McIV", "WingPh1", "WingPh2",
                             "WingPh3", "WingPh4", "Femur", "Tibia", "MtIII"),
                           colnames(POSTCRANIAL))

for (el in PROP_ELEMENTS) {
  POSTCRANIAL[[paste0("prop_", el)]] <- POSTCRANIAL[[el]] / POSTCRANIAL$SkeletalSize
}
PROP_VARS <- paste0("prop_", PROP_ELEMENTS)

prop_res <- do.call(rbind, lapply(TARGETS, function(tg) {
  do.call(rbind, lapply(PROP_VARS, function(pv) {
    x <- POSTCRANIAL[[pv]]; y <- POSTCRANIAL[[tg]]
    ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
    if (sum(ok) < MIN_FIT) return(NULL)
    m <- lm(log(y[ok]) ~ log(x[ok]))
    cf <- summary(m)$coefficients
    data.frame(target = tg, proportion = pv, n = sum(ok),
               slope = round(cf[2, 1], 4),
               r_squared = round(summary(m)$r.squared, 4),
               p_value = signif(cf[2, 4], 4), stringsAsFactors = FALSE)
  }))
}))

if (!is.null(prop_res)) {
  prop_res$p_adjusted_BH <- signif(p.adjust(prop_res$p_value, "BH"), 4)
  prop_res$significant <- prop_res$p_adjusted_BH < 0.05
  prop_res <- prop_res[order(-prop_res$r_squared), ]

  cat("Strongest proportion-metric relationships:\n")
  # Show a fixed number of decimals so every row lines up.
  prop_res_display <- prop_res
  prop_res_display$slope         <- format_fixed(prop_res$slope, 4)
  prop_res_display$r_squared     <- format_fixed(prop_res$r_squared, 4)
  prop_res_display$p_adjusted_BH <- format_fixed(prop_res$p_adjusted_BH, 4)
  print(head(prop_res_display[, c("target", "proportion", "n", "slope", "r_squared",
                          "p_adjusted_BH")], 12), row.names = FALSE)
  write_supp(prop_res_display, "S15_metric_vs_proportions")
}

cat("\nNOTE: independence is assumed throughout. With Pagel's lambda near 0.7,\n")
cat("refit any slope carried into the manuscript with a PGLS.\n\n")
