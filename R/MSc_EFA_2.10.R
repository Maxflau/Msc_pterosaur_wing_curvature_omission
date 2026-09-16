# ── 10H. SAVE ─────────────────────────────────────────────────────────────────
n_rows_efa <- ceiling(n_bins_efa / 3)
fig_h_efa  <- max(8, n_rows_efa * 5.2)

dir.create("clade_analyses_v4/plots",     showWarnings=FALSE, recursive=TRUE)
dir.create("clade_analyses_v4/plots_PDF", showWarnings=FALSE, recursive=TRUE)

ggsave("clade_analyses_v4/plots_PDF/CLADE_EFA_morphospace_time.pdf",
       efa_time_plot, width=18, height=fig_h_efa, dpi=300)
ggsave("clade_analyses_v4/plots/CLADE_EFA_morphospace_time.png",
       efa_time_plot, width=18, height=fig_h_efa, dpi=150)
cat("  ✓ EFA morphospace through time (PNG + PDF)\n")

cat("\n=== Section 11: EFA grid biomechanical table ===\n")

if (!exists("shape_perf")) {
  cat("  shape_perf not found - skipping.\n")
} else {
  
  if (!exists("res_dir")) res_dir <- "clade_analyses_v4/results"
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("  Specimens: %d\n", nrow(shape_perf)))
  
  # ── 11A. EXACT PC1/PC2 VIA ROTATION MATRIX ──────────────────────────────────
  PCA_VARS <- c("aspect_ratio","r2_hat","wing_loading",
                "von_mises_stress","wing_curvature","shape_complexity")
  
  missing_pca <- setdiff(PCA_VARS, names(shape_perf))
  if (length(missing_pca) > 0) {
    stop("Missing metrics: ", paste(missing_pca, collapse = ", "),
         "\n  Re-run Msc_performance_metrics_v3.R (with calculate_von_mises_stress sourced).")
  }
  
  # The rotation must have been fitted on the SAME six metrics. A loadings file
  # written by the old pipeline lists second_moment and stress_root, so
  # projecting the new metrics through it would silently produce wrong scores.
  use_loadings <- FALSE
  if (file.exists("data/pca_loadings.csv")) {
    loadings_df <- read.csv("data/pca_loadings.csv", stringsAsFactors = FALSE)
    rownames(loadings_df) <- loadings_df$variable
    if (all(PCA_VARS %in% rownames(loadings_df))) {
      use_loadings <- TRUE
    } else {
      cat("  pca_loadings.csv predates the metric rewrite - ignored.\n")
    }
  }
  
  if (use_loadings) {
    loadings_mat <- as.matrix(loadings_df[PCA_VARS, c("PC1","PC2")])
    ref_means <- sapply(PCA_VARS, function(v) mean(shape_perf[[v]], na.rm = TRUE))
    ref_sds   <- sapply(PCA_VARS, function(v) sd(shape_perf[[v]],   na.rm = TRUE))
    z_mat <- sweep(sweep(as.matrix(shape_perf[, PCA_VARS]), 2, ref_means, "-"),
                   2, ref_sds, "/")
    pc_proj  <- z_mat %*% loadings_mat
    pc1_vals <- round(pc_proj[, 1], 4)
    pc2_vals <- round(pc_proj[, 2], 4)
    cat("  PC1/PC2 from the rotation matrix.\n")
  } else if (all(c("PC1","PC2") %in% names(shape_perf))) {
    pc1_vals <- round(shape_perf$PC1, 4)
    pc2_vals <- round(shape_perf$PC2, 4)
    cat("  PC1/PC2 taken from shape_perf.\n")
  } else {
    stop("No usable loadings and no PC1/PC2 in shape_perf.")
  }
  
  
  join_missing <- function(target, path, cols, key = "species") {
    if (!file.exists(path)) return(target)
    src <- read.csv(path, stringsAsFactors = FALSE)
    if (!key %in% names(src)) return(target)
    src[[key]] <- trimws(as.character(src[[key]]))
    want <- intersect(setdiff(cols, names(target)), names(src))
    if (length(want) == 0) return(target)
    m <- match(trimws(as.character(target[[key]])), src[[key]])
    for (v in want) target[[v]] <- src[[v]][m]
    cat(sprintf("    from %s: %s (%d matched)\n", basename(path),
                paste(want, collapse = ", "), sum(!is.na(m))))
    target
  }
  
  cat("  Recovering missing columns:\n")
  
  # 1. Taxonomy and flight category, from the master database
  if (file.exists("data/pteros_main_data.csv")) {
    raw_db <- read.csv("data/pteros_main_data.csv", sep = ";", stringsAsFactors = FALSE)
    sp_col <- grep("^SPECIES", names(raw_db), value = TRUE)[1]
    if (!is.na(sp_col)) {
      raw_db$species <- trimws(raw_db[[sp_col]])
      m <- match(trimws(shape_perf$species), raw_db$species)
      for (v in c("family", "Flight.category", "Flight_category")) {
        tgt <- sub("Flight\\.category", "Flight_category", v)
        if (v %in% names(raw_db) && !tgt %in% names(shape_perf)) {
          shape_perf[[tgt]] <- trimws(as.character(raw_db[[v]][m]))
          cat(sprintf("    from pteros_main_data.csv: %s (%d matched)\n",
                      tgt, sum(!is.na(m))))
        }
      }
    }
  }
  
  # 2. Pareto and classification columns, from whichever CSV holds them
  PARETO_COLS <- c("goldberg_rank","strategy","limiting_objective",
                   "distance_to_front","on_front","pareto_rank_ratio")
  for (p in c("clade_analyses_v4/results/pareto_all_specimens.csv",
              "clade_analyses_v4/results/wing_classification.csv",
              "supplementals/wing_classification.csv",
              "clade_analyses_v4/results/reference_specimens_for_curvature.csv")) {
    shape_perf <- join_missing(shape_perf, p, PARETO_COLS)
  }
  
  # 3. Diet combination, derived rather than read
  if (!"Diet_combo" %in% names(shape_perf) &&
      all(c("Diet.1","Diet.2") %in% names(shape_perf))) {
    d1 <- trimws(as.character(shape_perf$Diet.1))
    d2 <- trimws(as.character(shape_perf$Diet.2))
    # Sorted pair, so "A+B" and "B+A" are one category: the raw columns contain
    # both orderings, which otherwise splits the counts across duplicates
    shape_perf$Diet_combo <- mapply(function(a, b) {
      if (is.na(a) || a == "") return(NA_character_)
      if (is.na(b) || b == "" || b == a) a else paste(sort(c(a, b)), collapse = "+")
    }, d1, d2, USE.NAMES = FALSE)
    cat(sprintf("    derived: Diet_combo (%d levels)\n",
                length(unique(na.omit(shape_perf$Diet_combo)))))
  }
  
  still_missing <- setdiff(c("family","Flight_category","goldberg_rank","strategy",
                             "limiting_objective","distance_to_front","Diet_combo"),
                           names(shape_perf))
  if (length(still_missing) > 0) {
    cat("  Still absent:", paste(still_missing, collapse = ", "), "\n")
    if (any(c("strategy","limiting_objective","distance_to_front") %in% still_missing)) {
      cat("    -> run Msc_wing_classification_v1.R\n")
    }
    if ("goldberg_rank" %in% still_missing) cat("    -> run Msc_pareto_front_v3.R\n")
  }
  # ── 11B. BUILD TABLE ────────────────────────────────────────────────────────
  ref_cols <- c("species","aspect_ratio","r2_hat","wing_loading",
                "von_mises_stress","wing_curvature","shape_complexity",
                "Order","clade","family","Environment","Diet.1","Diet.2",
                "Midpoint","Period.Name",
                "Depositional.settings.paleoenvironment","Flight_category",
                "pareto_rank_ratio","goldberg_rank","on_front","strategy",
                "limiting_objective","distance_to_front",
                "Diet_combo","Time_Bin")
  
  base_cols <- ref_cols[ref_cols %in% names(shape_perf)]
  
  # Report what is absent instead of dropping it silently: the previous version
  # wrote 18 of 34 requested columns without saying so
  dropped <- setdiff(ref_cols, base_cols)
  if (length(dropped) > 0) {
    cat("  Columns not available:", paste(dropped, collapse = ", "), "\n")
  }
  
  sp_table <- shape_perf[, base_cols, drop = FALSE]
  sp_table$shapePC1 <- round(shape_perf$shapePC1, 4)
  sp_table$shapePC2 <- round(shape_perf$shapePC2, 4)
  sp_table$PC1 <- pc1_vals
  sp_table$PC2 <- pc2_vals
  
  first_cols <- c("species","aspect_ratio","r2_hat","wing_loading",
                  "von_mises_stress","wing_curvature","shape_complexity",
                  "Order","clade","family","Environment","Diet.1","Diet.2",
                  "Midpoint","Period.Name",
                  "Depositional.settings.paleoenvironment","Flight_category",
                  "shapePC1","shapePC2","PC1","PC2")
  first_cols <- first_cols[first_cols %in% names(sp_table)]
  sp_table <- sp_table[, c(first_cols, setdiff(names(sp_table), first_cols))]
  
  # ── 11C. EXPORT ─────────────────────────────────────────────────────────────
  out_path <- file.path(res_dir, "specimen_biomechanical_pca_table.csv")
  write.csv(sp_table, out_path, row.names = FALSE)
  cat(sprintf("  %d specimens x %d columns -> %s\n",
              nrow(sp_table), ncol(sp_table), out_path))
  
  # Loadings written alongside, so the next run can project exactly
  if (exists("pca_performance")) {
    write.csv(data.frame(variable = rownames(pca_performance$rotation),
                         pca_performance$rotation[, 1:2]),
              "data/pca_loadings.csv", row.names = FALSE)
    cat("  pca_loadings.csv refreshed from the current PCA.\n")
  }
}