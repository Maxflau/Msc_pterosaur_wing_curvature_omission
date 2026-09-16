# Build the shared data objects that all later statistics scripts use.
cat("\n================================================================================\n")
cat("STATISTICS - SHARED SETUP\n")
cat("================================================================================\n\n")

library(dplyr)

OUT_DIR <- "output/results/supplementals"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Groups below this size are excluded: an SD or a variance on two points is
# meaningless, and a rank test on them is noise.
MIN_N  <- 3
N_PERM <- 9999
set.seed(2026)

# Save one table to the supplement folder and print what was written.
write_supp <- function(x, name) {
  if (is.null(x) || nrow(x) == 0) {
    cat(sprintf("  SKIPPED %s - no rows\n", name)); return(invisible(NULL))
  }
  write.csv(x, file.path(OUT_DIR, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("  Exported: %s/%s.csv (%d rows)\n", OUT_DIR, name, nrow(x)))
}

# --------------------------------------------------------------------------------
# Check that the shared table exists, then list the metrics and groups to use
# --------------------------------------------------------------------------------
if (!exists("performance_data_clean")) {
  stop("performance_data_clean not found - source the pipeline up to ", "Msc_impossible_region_v3.R first.")
}

#################################################################################
# Add the stress column here if an earlier script has not made it yet
#################################################################################
# If the stress column is missing, rebuild it here instead of skipping it.
# from the outlines. calculate_von_mises_stress() needs
# calculate_second_moment_DIAGNOSTIC() from the same functions script.
if (!"von_mises_stress" %in% colnames(performance_data_clean)) {
  if (exists("calculate_von_mises_stress") && exists("outlines_list")) {
    sv <- vapply(seq_len(nrow(performance_data_clean)), function(i) {
      s <- performance_data_clean$species[i]
      o <- outlines_list[[s]]
      if (is.null(o) || nrow(o) < 4) return(NA_real_)
      mk <- if ("Mass_kg" %in% colnames(performance_data_clean))
        performance_data_clean$Mass_kg[i] else NA_real_
      ws <- if ("Wingspan_cm" %in% colnames(performance_data_clean))
        performance_data_clean$Wingspan_cm[i] else NA_real_
      tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
    }, numeric(1))

    if (sum(is.finite(sv)) > 10) {
      performance_data_clean$von_mises_stress <- sv
      cat(sprintf("von_mises_stress computed here for %d / %d specimens\n",
                  sum(is.finite(sv)), length(sv)))
    } else {
      cat("von_mises_stress could not be computed - too few finite values.\n")
    }
  } else {
    missing_dep <- c(
      if (!exists("outlines_list")) "outlines_list [Msc_outline_loading.R]")
    cat("von_mises_stress missing, and cannot be computed. Absent:",
        paste(missing_dep, collapse = ", "), "\n")
  }
}

# The scale-invariant set. second_moment is deliberately absent; keep it out
# of every downstream test.
METRICS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress", "wing_curvature",
 "shape_complexity", "wing_loading_ratio", "pareto_rank_ratio"),colnames(performance_data_clean))

#################################################################################
# Build one combined diet label so every later script uses the same wording
#################################################################################

if (all(c("Diet.1", "Diet.2") %in% colnames(performance_data_clean))) {
  d1 <- trimws(as.character(performance_data_clean$Diet.1))
  d2 <- trimws(as.character(performance_data_clean$Diet.2))
  performance_data_clean$Diet_combined <- ifelse(
    is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = " + "))

  tb <- table(performance_data_clean$Diet_combined)
  cat(sprintf("\nDiet combinations: %d distinct, %d with n >= %d\n",
              length(tb), sum(tb >= MIN_N), MIN_N))
  print(sort(tb, decreasing = TRUE))
}

# Clade and Order are the two taxonomic scales; the rest are ecological
GROUPS  <- intersect(c("clade", "Order"), colnames(performance_data_clean))
FACTORS <- intersect(c("clade", "Order", "Diet.1", "Diet.2","Diet_combined",
                       "Depositional.settings.paleoenvironment", "Environment",
                       "Period.Name"),
                     colnames(performance_data_clean))

cat("Metrics:", paste(METRICS, collapse = ", "), "\n")
cat("Taxonomic levels:", paste(GROUPS, collapse = ", "), "\n")
cat("All factors:", paste(FACTORS, collapse = ", "), "\n")
cat("Specimens:", nrow(performance_data_clean), "| min group size:", MIN_N, "\n")

if (length(METRICS) < 2) {
  stop("Too few valid metrics - run Msc_performance_metrics_v3.R with von_mises_stress.")
}
if (!"von_mises_stress" %in% METRICS) {
  cat("\nWARNING: running WITHOUT a stress metric.\n")
  cat("To fix, source in this order before the statistics:\n")
  cat("  source(\"<biomechanical functions script>\")   # calculate_von_mises_stress()\n")
  cat("  source(\"Msc_outline_loading.R\")               # outlines_list\n")
  cat("Then re-run this setup: it will compute the column itself.\n")
}

#################################################################################
# Build the cleaned multivariate matrix used by the later group tests
#################################################################################
# Log then standardise, as for the performance PCA: the metrics have
# heterogeneous units, and a Euclidean distance on raw values would be dominated
# by wing loading (SD 26.8 against 0.11 for curvature).
MV_METRICS <- setdiff(METRICS, "pareto_rank_ratio")
Mraw <- as.matrix(performance_data_clean[, MV_METRICS])
mv_keep <- complete.cases(Mraw) & apply(Mraw > 0, 1, all)

MV_MATRIX <- scale(log(Mraw[mv_keep, , drop = FALSE]))
MV_META   <- performance_data_clean[mv_keep, ]

cat(sprintf("Multivariate matrix: %d specimens x %d metrics\n\n",
            nrow(MV_MATRIX), ncol(MV_MATRIX)))
cat("Setup complete. Available: METRICS, GROUPS, FACTORS, MV_MATRIX, MV_META,\n")
cat("write_supp(), MIN_N, N_PERM, OUT_DIR\n\n")
