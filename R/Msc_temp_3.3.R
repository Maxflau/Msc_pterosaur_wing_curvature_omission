FIG_W        <- 40     # width of a 3-column figure
FIG_H_ROW    <- 6    # height per row of panels
FIG_DPI      <- 600
PANEL_BORDER <- 0.9    # frame thickness around each panel


VAR_LABELS <- c(
  aspect_ratio      = "Aspect ratio",
  r2_hat            = "Second moment of area",
  von_mises_stress  = "Von Mises stress (scale-invariant)",
  wing_curvature    = "Wing curvature (camber)",
  shape_complexity  = "Shape complexity",
  pareto_score_aspect_ratio     = "Pareto score - Aspect ratio",
  pareto_score_r2_hat           = "Pareto score - Second moment of area",
  pareto_score_von_mises_stress = "Pareto score - Von Mises stress",
  Wingspan_cm       = "Wingspan (cm, log10)",
  Mass_kg           = "Body mass (kg, log10)",
  Wing_area_cm2     = "Wing area (cm2, log10)",
  pareto_rank_ratio = "Pareto rank ratio")

# One colour per trait, so no two panels share a hue
VAR_COLOURS <- c(
  aspect_ratio      = "#1F6FB4",   # blue
  r2_hat            = "#B5651D",   # ochre
  von_mises_stress  = "#8E3B8E",   # purple
  wing_curvature    = "#2E8B57",   # green
  shape_complexity  = "#C1443C",   # red
  pareto_score_aspect_ratio     = "#1F6FB4",
  pareto_score_r2_hat           = "#B5651D",
  pareto_score_von_mises_stress = "#8E3B8E",
  Wingspan_cm       = "#B5651D",
  Mass_kg           = "#8E3B8E",
  Wing_area_cm2     = "#2E8B57",
  pareto_rank_ratio = "#1b7837")

CAPTION <- paste("Mean with bootstrap 95% percentile interval (2000 resamples).",
                 "Specimens within a bin are not phylogenetically independent,",
                 "so intervals are optimistic.")

# --------------------------------------------------------------------------------
# Panel builder
# --------------------------------------------------------------------------------
make_time_panel <- function(vars, title, subtitle, ncol = 3) {
  
  d <- temporal_metrics[temporal_metrics$variable %in% vars, ]
  if (nrow(d) == 0) return(NULL)
  
  keep_vars <- vars[vars %in% d$variable]
  d$label <- factor(VAR_LABELS[d$variable], levels = VAR_LABELS[keep_vars])
  d$col   <- VAR_COLOURS[d$variable]
  
  ggplot(d, aes(x = Time_Bin, y = mean, group = 1)) +
    # Shaded band and error bar show the same interval two ways: the band makes
    # the trend readable, the bar keeps each bin's uncertainty explicit
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = label), alpha = 0.18,
                colour = NA) +
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = label),
                  width = 0.16, linewidth = 0.55) +
    geom_line(aes(colour = label), linewidth = 0.85) +
    geom_point(aes(colour = label), size = 2.7) +
    geom_text(aes(label = paste0("n=", n)), y = -Inf, vjust = -0.7,
              size = 2.5, colour = "grey35") +
    scale_colour_manual(values = setNames(VAR_COLOURS[keep_vars],
                                          VAR_LABELS[keep_vars]), guide = "none") +
    scale_fill_manual(values = setNames(VAR_COLOURS[keep_vars],
                                        VAR_LABELS[keep_vars]), guide = "none") +
    facet_wrap(~ label, scales = "free_y", ncol = ncol) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey95", colour = "grey40"),
          panel.border = element_rect(colour = "grey30", fill = NA,
                                      linewidth = PANEL_BORDER),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8.5, colour = "grey35"))
}

save_fig <- function(p, name, n_panels, ncol = 3) {
  if (is.null(p)) return(invisible(NULL))
  h <- FIG_H_ROW * ceiling(n_panels / ncol)
  ggsave(sprintf("clade_analyses_v4/plots/%s.png", name), p,
         width = FIG_W, height = h, dpi = FIG_DPI)
  ggsave(sprintf("clade_analyses_v4/plots_PDF/%s.pdf", name), p,
         width = FIG_W, height = h)
  cat(sprintf("Written: %s (%.0f x %.1f in, %d panels)\n", name, FIG_W, h, n_panels))
}

present <- unique(temporal_metrics$variable)
stress_var <- intersect(c("von_mises_stress"), present)[1]

# --------------------------------------------------------------------------------
# Failsafe: build the stress rows here if they are absent
# --------------------------------------------------------------------------------
# The stress panel kept vanishing because von_mises_stress was missing from
# temporal_metrics, and intersect() dropped it silently. Rather than depend on
# the upstream chain, compute the bin means directly from the outlines.
if (is.na(stress_var) && exists("temporal_data") && exists("outlines_list") &&
    exists("calculate_von_mises_stress")) {
  
  sv <- vapply(seq_len(nrow(temporal_data)), function(i) {
    s <- temporal_data$species[i]
    o <- outlines_list[[s]]
    if (is.null(o)) return(NA_real_)
    mk <- if ("Mass_kg" %in% colnames(temporal_data)) temporal_data$Mass_kg[i] else NA_real_
    ws <- if ("Wingspan_cm" %in% colnames(temporal_data)) temporal_data$Wingspan_cm[i] else NA_real_
    tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
  }, numeric(1))
  
  if (sum(is.finite(sv)) > 10) {
    boot_ci <- function(x, n_boot = 2000) {
      x <- x[is.finite(x)]
      if (length(x) < 3) return(c(NA, NA, NA, length(x)))
      bm <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
      c(mean(x), unname(quantile(bm, 0.025)), unname(quantile(bm, 0.975)), length(x))
    }
    add <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
      ci <- boot_ci(sv[temporal_data$Time_Bin == tb])
      data.frame(Time_Bin = tb, variable = "von_mises_stress", scale = "raw",
                 mean = ci[1], lo = ci[2], hi = ci[3], n = ci[4],
                 stringsAsFactors = FALSE)
    }))
    add$Time_Bin <- factor(add$Time_Bin, levels = levels(temporal_metrics$Time_Bin))
    temporal_metrics <- rbind(temporal_metrics, add)
    present <- unique(temporal_metrics$variable)
    stress_var <- "von_mises_stress"
    cat(sprintf("von_mises_stress bin means computed here for %d specimens\n",
                sum(is.finite(sv))))
  }
}

###################################################################################
# 1. Pareto score of the three Pareto objectives (aspect ratio, second moment
#    of area, von Mises stress) - each metric's within-sample percentile
#    (0-1), not its raw value. From Msc_pareto_front_1.R.
###################################################################################
shape_vars <- intersect(c("pareto_score_aspect_ratio", "pareto_score_r2_hat","pareto_score_von_mises_stress"), present)

cat("Shape panel variables:", paste(shape_vars, collapse = ", "), "\n")
if (length(shape_vars) == 0) {
  cat("NOTE: pareto_score_* columns absent from temporal_metrics. Source\n")
  cat("Msc_pareto_front_1.R before Msc_temp_3_1.R so they get bootstrapped.\n")
}
p_shape <- make_time_panel(shape_vars, "Pareto score of the Pareto objectives through the Mesozoic", CAPTION)
save_fig(p_shape, "TEMPORAL_01_shape_traits", length(shape_vars))