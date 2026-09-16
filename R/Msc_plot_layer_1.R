##################################################################################
# 0. Constants, guarded
##################################################################################
if (!exists("performance_data_clean")) {
  stop("performance_data_clean not found - source Msc_impossible_region_v3.R first.")
}
if (!exists("inside_envelope")) {
  stop("inside_envelope() not found - source Msc_impossible_region_v3.R first.")
}

# ── 1. ENSURE von_mises_stress EXISTS ─────────────────────────────────────────
if (!"von_mises_stress" %in% colnames(performance_data_clean)) {
  
  cat("von_mises_stress absent - attempting to build it.\n")
  
  # (a) from the outlines, if the function and the outlines are in memory
  if (exists("calculate_von_mises_stress") && exists("outlines_list")) {
    performance_data_clean$von_mises_stress <- vapply(
      seq_len(nrow(performance_data_clean)), function(i) {
        s <- performance_data_clean$species[i]
        o <- outlines_list[[s]]
        if (is.null(o) || nrow(o) < 4) return(NA_real_)
        mk <- if ("Mass_kg" %in% colnames(performance_data_clean))
          performance_data_clean$Mass_kg[i] else NA_real_
        ws <- if ("Wingspan_cm" %in% colnames(performance_data_clean))
          performance_data_clean$Wingspan_cm[i] else NA_real_
        tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
      }, numeric(1))
    cat(sprintf("  computed from outlines: %d of %d specimens\n",
                sum(is.finite(performance_data_clean$von_mises_stress)),
                nrow(performance_data_clean)))
  }
  
  # (b) from a CSV written by an earlier run
  if (!"von_mises_stress" %in% colnames(performance_data_clean) ||
      all(is.na(performance_data_clean$von_mises_stress))) {
    for (p in c("output/results/pareto_all_specimens.csv",
                "output/results/performance_metrics_complete_CLADE.csv",
                "output/results/supplementals/wing_classification.csv")) {
      if (!file.exists(p)) next
      src <- read.csv(p, stringsAsFactors = FALSE)
      if (!all(c("species", "von_mises_stress") %in% names(src))) next
      m <- match(trimws(performance_data_clean$species), trimws(src$species))
      performance_data_clean$von_mises_stress <- src$von_mises_stress[m]
      cat(sprintf("  read from %s: %d matched\n", basename(p), sum(!is.na(m))))
      break
    }
  }
  
  # (c) give up loudly rather than silently plotting another metric
  if (!"von_mises_stress" %in% colnames(performance_data_clean) ||
      sum(is.finite(performance_data_clean$von_mises_stress)) < 20) {
    stop("von_mises_stress could not be built.\n",
         "  Source in this order, then re-run this file:\n",
         "    Msc_biomechanical_functions_v3.R  (calculate_von_mises_stress)\n",
         "    Msc_outline_loading.R             (outlines_list)\n",
         "    Msc_performance_metrics_v3.R")
  }
}

# ── 2. BACKGROUND VARIABLE ────────────────────────────────────────────────────
if (!exists("BACKGROUND_VAR")) BACKGROUND_VAR <- "von_mises_stress"
stopifnot(is.character(BACKGROUND_VAR), length(BACKGROUND_VAR) == 1)

# No silent fallback: a surface showing a different metric from the one named in
# the title is worse than an error
if (!BACKGROUND_VAR %in% colnames(performance_data_clean)) {
  stop("BACKGROUND_VAR '", BACKGROUND_VAR, "' not in performance_data_clean. ",
       "Available: ", paste(intersect(c("von_mises_stress","r2_hat","aspect_ratio",
                                        "wing_loading", "wing_curvature","shape_complexity","pareto_rank_ratio"),
                                      colnames(performance_data_clean)), collapse = ", "))
}

BACKGROUND_LABELS <- c(
  von_mises_stress  = "Von Mises stress\n(scale-invariant)",
  r2_hat            = "Second moment\nof area",
  aspect_ratio      = "Aspect ratio",
  wing_curvature    = "Wing curvature",
  shape_complexity  = "Shape complexity",
  wing_loading      = "Wing loading\n(N/m2)",
  pareto_rank_ratio = "Pareto\nrank ratio")

surf_label <- if (BACKGROUND_VAR %in% names(BACKGROUND_LABELS)) {
  BACKGROUND_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR

cat(sprintf("\nSurface variable: %s | colour bar: %s\n", BACKGROUND_VAR, gsub("\n", " ", surf_label)))

#######################################################################
#BACKGROUND GRID
######################################################################

pc1_rng <- range(performance_data_clean$PC1, na.rm = TRUE)
pc2_rng <- range(performance_data_clean$PC2, na.rm = TRUE)
pad1 <- diff(pc1_rng) * 0.15; pad2 <- diff(pc2_rng) * 0.15

stress_grid <- expand.grid( PC1 = seq(pc1_rng[1] - pad1, pc1_rng[2] + pad1, length.out = 100), PC2 = seq(pc2_rng[1] - pad2, pc2_rng[2] + pad2, length.out = 100))

# Bandwidth scaled to the morphospace, not a hard-coded 0.5 that changes meaning
# whenever the PCA changes
bw <- 0.08 * max(diff(pc1_rng), diff(pc2_rng))
zvals <- performance_data_clean[[BACKGROUND_VAR]]

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

# Cells far beyond the padded envelope (bw's Gaussian kernel underflows to
# zero there) fall back to the single nearest specimen's value rather than
# NA, so the colour fades out smoothly toward the edge of the grid instead
# of breaking to a blank white patch.
stress_grid$stress <- vapply(seq_len(nrow(stress_grid)), function(i) {
  d2 <- (performance_data_clean$PC1 - stress_grid$PC1[i])^2 +
    (performance_data_clean$PC2 - stress_grid$PC2[i])^2
  w <- exp(-d2 / (2 * bw^2))
  if (sum(w) > 1e-6) {
    weighted_median(zvals, w)
  } else {
    zvals[which.min(d2)]
  }
}, numeric(1))

stress_grid$viable <- inside_envelope(stress_grid$PC1, stress_grid$PC2)

cat(sprintf("Grid: %d cells | %d inside the envelope\n",
            nrow(stress_grid), sum(stress_grid$viable)))
cat(sprintf("Kernel bandwidth: %.4g | weighted median, nearest-neighbour fallback\n", bw))

######################################################################
#COLOUR SCALE 
######################################################################

sv <- stress_grid$stress[is.finite(stress_grid$stress)]
if (length(sv) == 0) stop("Background grid is entirely NA.")

surf_limits <- unname(quantile(sv, c(0.01, 0.99)))
if (surf_limits[1] == surf_limits[2]) surf_limits <- range(sv)
use_log <- all(sv > 0) && (surf_limits[2] / surf_limits[1] > 10)

STRESS_PALETTES <- list(
  dark  = c("#6baed6","#4292c6","#2171b5","#08519c","#08306b","#03152e"),
  light = c("#f7fbff","#deebf7","#c6dbef","#9ecae1", "#6baed6","#4292c6","#2171b5","#08519c"))

#' Drop-in replacement used by all the plot scripts.
scale_fill_stress <- function(palette = "dark", name = surf_label) {
  palette <- match.arg(palette, names(STRESS_PALETTES))
  ggplot2::scale_fill_gradientn( colors   = STRESS_PALETTES[[palette]],
                                 name     = name, na.value = "white",
                                 trans    = if (use_log) "log10" else "identity",
                                 limits   = surf_limits,
                                 oob      = scales::squish,
                                 guide    = ggplot2::guide_colourbar(
                                   barheight      = grid::unit(60, "mm"), barwidth = grid::unit(5,  "mm"),
                                   ticks.colour   = "grey30",frame.colour   = "grey30", title.position = "top",
                                   label.theme    = ggplot2::element_text(size = 9, hjust = 0)))
}

cat(sprintf("Scale: %s | limits %.4g to %.4g\n\n",
            if (use_log) "log10" else "linear", surf_limits[1], surf_limits[2]))

ENVELOPE_CAPTION <- paste0("Surface: ", gsub("\n", " ", surf_label), ".")