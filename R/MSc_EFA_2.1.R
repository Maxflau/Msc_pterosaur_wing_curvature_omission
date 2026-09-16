select <- dplyr::select; filter <- dplyr::filter
CLADE_COLS <- c(
  "Anurognathidae"         = "#46852F",
  "Azhdarchoidea"          = "#DB8E00",
  "Basal pterodactyliform" = "#F4A7B9",
  "Basal pterosaur"        = "#0CE27F",
  "Ctenochasmatoidea"      = "#9F1212",
  "Darwinoptera"           = "#791444",
  "Dsungaripteroidea"      = "#60A293",
  "Eudimorphodontoidea"    = "#FF00F0",
  "Ornithocheiromorpha"    = "#FFE413",
  "Pteranodontia"          = "#DA8AAA",
  "Rhamphorhynchidae"      = "#03240B")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
# Load the main results table and rebuild key fields if they are missing.
PERF_CSV <- "output/results/performance_metrics_complete_CLADE.csv"

if (exists("performance_data_clean")) {
  perf <- performance_data_clean
  cat(sprintf("perf: using performance_data_clean in memory (%d rows)\n", nrow(perf)))
} else if (file.exists(PERF_CSV)) {
  perf <- read.csv(PERF_CSV, stringsAsFactors = FALSE)
  cat(sprintf("perf: read from CSV (%d rows) - may predate the metric rewrite\n",
              nrow(perf)))
} else {
  stop("Neither performance_data_clean nor ", PERF_CSV, " is available.")
}

# Recompute von_mises_stress when the source predates it
if (!"von_mises_stress" %in% names(perf)) {
  if (exists("calculate_von_mises_stress") && exists("outlines_list")) {
    perf$von_mises_stress <- vapply(seq_len(nrow(perf)), function(i) {
      s <- perf$species[i]
      o <- outlines_list[[s]]
      if (is.null(o) || nrow(o) < 4) return(NA_real_)
      mk <- if ("Mass_kg" %in% colnames(perf)) perf$Mass_kg[i] else NA_real_
      ws <- if ("Wingspan_cm" %in% colnames(perf)) perf$Wingspan_cm[i] else NA_real_
      tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
    }, numeric(1))
    cat(sprintf("von_mises_stress recomputed for %d specimens\n",
                sum(is.finite(perf$von_mises_stress))))
  } else {
    stop("von_mises_stress absent and cannot be recomputed. Source in this order:\n",
         "  <biomechanical functions script> (calculate_von_mises_stress) ->\n",
         "  Msc_outline_loading.R -> Msc_performance_metrics_v3.R")
  }
}

for (v in c("r2_hat", "pareto_rank_ratio")) {
  if (!v %in% names(perf)) {
    stop(v, " absent from perf. Re-run ",
         ifelse(v == "r2_hat", "Msc_performance_metrics_v3.R",
                "Msc_pareto_front_v3.R"), " before this script.")
  }
}

perm_csv <- read.csv("output/results/04_permanova_pc1pc2.csv",
                     stringsAsFactors = FALSE)

get_perm <- function(factor_name) {
  r <- perm_csv[perm_csv$Factor == factor_name, ]
  if (nrow(r) == 0) return(list(r2 = NA, f = NA, p = NA, sig = ""))
  list(r2 = round(r$R2[1], 3), f = round(r$F[1], 1),
       p = r$p_value[1], sig = r$sig[1])
}

PERF_JOIN_COLS <- c("species","aspect_ratio","wing_loading","von_mises_stress",
                    "wing_curvature","shape_complexity","r2_hat",
                    "pareto_rank_ratio","on_front","strategy",
                    "PC1","PC2","clade","Order","Diet.1","Diet.2",
                    "Environment","Depositional.settings.paleoenvironment",
                    "Midpoint","Period.Name")

build_shape_perf <- function() {
  cols_to_add <- PERF_JOIN_COLS[PERF_JOIN_COLS %in% names(perf)]
  if (exists("shape_scores")) {
    cols_to_add <- c("species", cols_to_add[!cols_to_add %in% names(shape_scores)])
  }
  sp <- shape_scores %>%
    left_join(perf %>% select(any_of(cols_to_add)), by = "species") %>%
    select(-ends_with(".y")) %>%
    rename_with(~ str_remove(.x, "\\.x$")) %>%
    filter(!is.na(pareto_rank_ratio))

  cat(sprintf("shape_perf: %d specimens | on Pareto front: %d | von_mises_stress: %d\n",
              nrow(sp), sum(sp$on_front, na.rm = TRUE),
              sum(is.finite(sp$von_mises_stress))))
  sp
}

if (exists("shape_scores")) {
  shape_perf <- build_shape_perf()
} else {
  message("shape_scores not found - run build_shape_perf() after MSc_rep.R")
}
if (all(c("Diet.1", "Diet.2") %in% names(shape_perf)) &&
    !"Diet_combined" %in% names(shape_perf)) {
  d1 <- trimws(as.character(shape_perf$Diet.1))
  d2 <- trimws(as.character(shape_perf$Diet.2))
  shape_perf$Diet_combined <- ifelse(is.na(d2) | d2 == "" | d2 == d1,
                                     d1, paste(d1, d2, sep = "+"))
  cat(sprintf("Diet_combined: %d distinct levels, %d with n >= 3\n",
              length(unique(shape_perf$Diet_combined)),
              sum(table(shape_perf$Diet_combined) >= 3)))
  print(sort(table(shape_perf$Diet_combined), decreasing = TRUE))
}
# ── 3. TIME BINS ──────────────────────────────────────────────────────────────
# Make the time groups used by the later time-based plots.
# Section 10 skipped entirely because Time_Bin was absent from shape_perf.
if (!"Time_Bin" %in% names(shape_perf)) {
  pcol <- intersect(c("Period.Name", "Period_Name"), names(shape_perf))[1]
  if (!is.na(pcol)) {
    pc <- trimws(as.character(shape_perf[[pcol]]))
    pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
    shape_perf$Time_Bin <- factor(pc,
                                  levels = c("Late Triassic", "Early+Middle Jurassic", "Late Jurassic",
                                             "Early Cretaceous", "Late Cretaceous"))
    cat(sprintf("Time_Bin built: %d specimens across %d bins\n",
                sum(!is.na(shape_perf$Time_Bin)),
                nlevels(droplevels(shape_perf$Time_Bin))))
  } else {
    cat("No period column - Time_Bin not built, Section 10 will skip.\n")
  }
}
if (!"Midpoint" %in% names(shape_perf)) {
  cat("WARNING: Midpoint absent - Section 10 needs it as well as Time_Bin.\n")
}

# ── 4. DISTANCE MATRICES ──────────────────────────────────────────────────────
# Build the distance tables used by the later comparison tests.
# shape_dist existed; bio_dist did not, so betadisper failed on the second call.
if (!exists("shape_dist")) {
  shape_dist <- dist(as.matrix(shape_perf[, c("shapePC1", "shapePC2")]),
                     method = "euclidean")
}

if (!exists("bio_dist")) {
  BIO_FOR_DIST <- intersect(c("aspect_ratio", "r2_hat", "wing_loading",
                              "von_mises_stress", "wing_curvature", "shape_complexity"),
                            names(shape_perf))
  M <- as.matrix(shape_perf[, BIO_FOR_DIST])

  # Log then standardise, as for the performance PCA: without it the Euclidean
  # distance is dominated by wing loading (SD 26.8 against 0.11 for curvature)
  ok <- complete.cases(M) & apply(M > 0, 1, all)
  if (sum(ok) < nrow(M)) {
    cat(sprintf("bio_dist: %d of %d specimens have all metrics finite and positive.\n",
                sum(ok), nrow(M)))
    cat("Rows with missing metrics are excluded; betadisper needs a matching subset.\n")
  }
  bio_mat  <- scale(log(M[ok, , drop = FALSE]))
  bio_dist <- dist(bio_mat, method = "euclidean")
  bio_keep <- ok          # index to subset stat_df when using bio_dist
  cat(sprintf("bio_dist: %d specimens x %d metrics\n",
              sum(ok), length(BIO_FOR_DIST)))
}
