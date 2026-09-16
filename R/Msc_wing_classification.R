# Classify each taxon by which performance trade-off appears to hold it back.
OBJECTIVES <- c(aspect_ratio = "max", aspect_ratio = "min",   r2_hat = "min", r2_hat = "max", 
                von_mises_stress = "min", von_mises_stress = "max", wing_loading_ratio ="max", wing_loading_ratio = "min")
OBJ_MEANING <- c(
  aspect_ratio     = "aerodynamic efficiency (lift-to-drag)",
  r2_hat           = "distal area distribution (average lift in flapping)",
  von_mises_stress = "root bending load (scale-invariant)",
  wing_loading_ratio = "gliding speed")

# --------------------------------------------------------------------------------
# 0. Check which performance objectives can be used for the ranking
# --------------------------------------------------------------------------------
missing_obj <- setdiff(names(OBJECTIVES), colnames(performance_data_clean))

if (length(missing_obj) > 0) {
  cat("Missing objective column(s):", paste(missing_obj, collapse = ", "), "\n")
  cat("Present metrics:", paste(intersect(
    c("aspect_ratio", "r2_hat", "von_mises_stress",
      "wing_loading_ratio", "wing_curvature", "shape_complexity"),
    colnames(performance_data_clean)), collapse = ", "), "\n")

  OBJECTIVES <- OBJECTIVES[!names(OBJECTIVES) %in% missing_obj]
  cat("Dropped:", paste(missing_obj, collapse = ", "),
      "- running with", length(OBJECTIVES), "objectives.\n")

  if (length(OBJECTIVES) < 2) {
    stop("Fewer than two objectives available - a Pareto front needs at least ",
         "two. Source the biomechanical functions script (calculate_von_mises_stress) ",
         "and Msc_performance_metrics_v3.R first.")
  }
  cat("\n")
}

cat("Objectives in use:", paste(names(OBJECTIVES), OBJECTIVES,
                                sep = " (", collapse = "), "), ")\n\n")

# --------------------------------------------------------------------------------
# LEVEL 1 and 2 — sort theoretical shapes into broad outcome classes
# --------------------------------------------------------------------------------
if (exists("theoretical_data") &&
    all(names(OBJECTIVES) %in% colnames(theoretical_data))) {
  theoretical_data$geometry_class <- ifelse(
    theoretical_data$self_intersecting, "Geometrically impossible", "Possible")

  obj_ranges <- lapply(names(OBJECTIVES),
                       function(v) range(performance_data_clean[[v]], na.rm = TRUE))
  names(obj_ranges) <- names(OBJECTIVES)

  within_observed <- Reduce(`&`, lapply(names(OBJECTIVES), function(v) {
    x <- theoretical_data[[v]]
    is.finite(x) & x >= obj_ranges[[v]][1] & x <= obj_ranges[[v]][2]
  }))

  theoretical_data$realisation_class <- ifelse(
    theoretical_data$self_intersecting, "Geometrically impossible",
    ifelse(within_observed, "Realised region", "Possible but unrealised"))

  cat("THEORETICAL SHAPES:\n")
  print(table(theoretical_data$realisation_class))
  cat("\n'Possible but unrealised' is the non-circular result: shapes that could\n")
  cat("exist and that no pterosaur occupied. The cause - mechanical, phylogenetic\n")
  cat("or ecological - is NOT determined by this analysis.\n\n")
} else {
  cat("theoretical_data unavailable or missing objectives - levels 1-2 skipped.\n\n")
}

# --------------------------------------------------------------------------------
# LEVEL 3 — classify real taxa by the trade-off that limits them most
# --------------------------------------------------------------------------------
obj_matrix <- as.matrix(performance_data_clean[, names(OBJECTIVES), drop = FALSE])
keep <- complete.cases(obj_matrix) & apply(is.finite(obj_matrix), 1, all)

if (sum(keep) < 10) {
  stop("Only ", sum(keep), " specimens have complete objective data.")
}

cls_data <- performance_data_clean[keep, ]
M <- obj_matrix[keep, , drop = FALSE]

for (v in names(OBJECTIVES)) if (OBJECTIVES[[v]] == "min") M[, v] <- -M[, v]

# Count how many real specimens can be ranked with the chosen objectives.
cat(sprintf("Real taxa classified: %d\n\n", nrow(M)))

cls_data$dominance_count <- vapply(seq_len(nrow(M)), function(i) {
  ge <- sweep(M, 2, M[i, ], ">="); gt <- sweep(M, 2, M[i, ], ">")
  sum(apply(ge, 1, all) & apply(gt, 1, any))
}, numeric(1))
cls_data$on_front <- cls_data$dominance_count == 0

pct <- apply(M, 2, function(x) rank(x, ties.method = "average") / length(x))
colnames(pct) <- names(OBJECTIVES)

cls_data$limiting_objective  <- colnames(pct)[apply(pct, 1, which.min)]
cls_data$limiting_percentile <- apply(pct, 1, min)
cls_data$balance             <- apply(pct, 1, function(p) 1 - (max(p) - min(p)))

Ms <- scale(M)
front_s <- Ms[cls_data$on_front, , drop = FALSE]
cls_data$distance_to_front <- vapply(seq_len(nrow(Ms)), function(i) {
  min(sqrt(rowSums((sweep(front_s, 2, Ms[i, ], "-"))^2)))
}, numeric(1))

cls_data$strategy <- with(cls_data, ifelse(
  on_front, "Best trade-off",
  ifelse(balance > 0.65, "Middling on all axes",
         paste0("Held back by ", limiting_objective))))

cat("STRATEGY CLASSIFICATION:\n")
print(table(cls_data$strategy))

cat("\nLimiting objective across all taxa:\n")
print(table(cls_data$limiting_objective))
cat("\nMeaning of each objective:\n")
for (v in names(OBJECTIVES)) {
  if (v %in% names(OBJ_MEANING)) cat(sprintf("  %-14s %s\n", v, OBJ_MEANING[[v]]))
}
cat("\n")

for (grp in intersect(c("clade", "Diet.1"), colnames(cls_data))) {
  cat(sprintf("Limiting objective by %s:\n", grp))
  tb <- table(cls_data[[grp]], cls_data$limiting_objective)
  print(tb)
  if (nrow(tb) > 1 && ncol(tb) > 1) {
    ft <- fisher.test(tb, simulate.p.value = TRUE, B = 10000)
    cat(sprintf("Fisher exact (simulated): p = %s\n\n",
                formatC(ft$p.value, format = "e", digits = 3)))
  } else cat("\n")
}

if ("Wingspan_cm" %in% colnames(cls_data)) {
  cat("Wingspan by limiting objective (cm):\n")
  print(tapply(cls_data$Wingspan_cm, cls_data$limiting_objective,
               function(x) round(summary(x[is.finite(x)]), 1)))
  cat("\n")
}

cat("NOTE: these tests treat specimens as independent. With Pagel's lambda near\n")
cat("0.7 they are not - confirm with the phylANOVA before interpreting.\n\n")

new_cols <- c("dominance_count", "on_front", "limiting_objective",
              "limiting_percentile", "balance", "distance_to_front", "strategy")

for (nm in new_cols) {
  performance_data_clean[[nm]] <- NA
  performance_data_clean[[nm]][keep] <- cls_data[[nm]]
}

dir.create("output/results", showWarnings = FALSE, recursive = TRUE)
write.csv(cls_data[, intersect(c("species", "clade", "Diet.1", "Diet.2", "DIETCOMBO_COLS", "Wingspan_cm",
                                 names(OBJECTIVES), new_cols), colnames(cls_data))],
          "output/results/wing_classification.csv", row.names = FALSE)

cat("Added to performance_data_clean:", paste(new_cols, collapse = ", "), "\n")
cat("Written: output/results/wing_classification.csv\n\n")
