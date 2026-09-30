# Draw the first set of temporal summary plots from the prepared time-bin tables.
###################################################################################
# 1. Give each plotted variable a readable label
###################################################################################
if (!exists("VAR_LABELS")) VAR_LABELS <- c(
  aspect_ratio      = "Aspect ratio",
  r2_hat            = "Second moment of area (r2-hat)",
  von_mises_stress  = "Von Mises stress (scale-invariant)",
  stress_index      = "Bending stress index (beam)",
  wing_curvature    = "Wing curvature (camber)",
  shape_complexity  = "Shape complexity",
  pareto_score_aspect_ratio     = "Pareto score - Aspect ratio",
  pareto_score_r2_hat           = "Pareto score - Second moment of area",
  pareto_score_von_mises_stress = "Pareto score - Von Mises stress",
  Wingspan_cm       = "Wingspan (cm, log10)",
  Mass_kg           = "Body mass (kg, log10)",
  Wing_area_cm2     = "Wing area (cm2, log10)",
  pareto_rank_ratio = "Pareto rank ratio")

plot_df <- temporal_metrics
plot_df$label <- factor(VAR_LABELS[plot_df$variable], levels = VAR_LABELS[VAR_LABELS %in% VAR_LABELS[plot_df$variable]])

###################################################################################
# 2. Make one reusable plotting function for any chosen variable set
###################################################################################
make_time_panel <- function(vars, title, subtitle, ncol = 3) {
  d <- plot_df[plot_df$variable %in% vars, ]
  if (nrow(d) == 0) return(NULL)

  ggplot(d, aes(x = Time_Bin, y = mean, group = 1)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4292c6", alpha = 0.20) + geom_line(colour = "#08519c", linewidth = 0.9) +
    geom_point(colour = "#08519c", size = 2.6) + geom_text(aes(label = paste0("n=", n)), y = -Inf, vjust = -0.8, size = 2.7, colour = "grey40") +
    facet_wrap(~ label, scales = "free_y", ncol = ncol) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
                                          strip.text = element_text(face = "bold", size = 9), panel.grid.minor = element_blank())
}

CAPTION_CI <- paste("Bootstrap 95% percentile intervals.",
                    "Specimens within a bin are not phylogenetically independent,",
                    "so intervals are optimistic.")

###################################################################################
# 3. Plot the three main Pareto score traits through time
#    moment of area and von Mises stress - the three Pareto objectives.
###################################################################################
p_shape <- make_time_panel(
  intersect(c("pareto_score_aspect_ratio", "pareto_score_r2_hat",
              "pareto_score_von_mises_stress"),
            unique(plot_df$variable)),
  "Pareto score of the Pareto objectives through the Mesozoic",
  paste("Each metric's within-sample percentile (0 = worst, 1 = best),",
        "not its raw value.", CAPTION_CI))

if (!is.null(p_shape)) {
  ggsave("output/plots_PDF/TEMPORAL_01_shape_traits.pdf", p_shape, width = 11, height = 4)
  cat("Written: TEMPORAL_01_shape_traits\n")
}

###################################################################################
# 4. Plot the body size variables through time
###################################################################################
# Size is now its own panel rather than an unacknowledged driver of the shape
# metrics. If size rises across the Cretaceous while shape traits stay flat, that
# separation is itself the result.
p_size <- make_time_panel(
  intersect(c("Wingspan_cm", "Mass_kg", "Wing_area_cm2"), unique(plot_df$variable)),
  "Body size through the Mesozoic",
  paste("Log\u2081\u2080 scale.", CAPTION_CI))

if (!is.null(p_size)) {
  ggsave("output/plots_PDF/TEMPORAL_02_body_size.pdf", p_size, width = 11, height = 4)
  cat("Written: TEMPORAL_02_body_size\n")
}

# --------------------------------------------------------------------------------
# 5. Plot disparity and Pareto rank ratio through time
# --------------------------------------------------------------------------------
disp_long <- rbind(
  data.frame(Time_Bin = disparity_metrics$Time_Bin, n = disparity_metrics$n,
             label = "Shape disparity (sum of variance)",
             mean = disparity_metrics$sov,
             lo = disparity_metrics$sov_lo, hi = disparity_metrics$sov_hi),
  data.frame(Time_Bin = disparity_metrics$Time_Bin, n = disparity_metrics$n,
             label = "Morphospace expansion (mean pairwise distance)",
             mean = disparity_metrics$mpd,
             lo = disparity_metrics$mpd_lo, hi = disparity_metrics$mpd_hi))

if ("pareto_rank_ratio" %in% plot_df$variable) {
  prr <- plot_df[plot_df$variable == "pareto_rank_ratio", ]
  disp_long <- rbind(disp_long,
                     data.frame(Time_Bin = prr$Time_Bin, n = prr$n,
                                label = "Pareto rank ratio",
                                mean = prr$mean, lo = prr$lo, hi = prr$hi))
}

p_disp <- ggplot(disp_long, aes(x = Time_Bin, y = mean, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#5aae61", alpha = 0.20) +
  geom_line(colour = "#1b7837", linewidth = 0.9) +
  geom_point(colour = "#1b7837", size = 2.6) +
  geom_text(aes(label = paste0("n=", n)), y = -Inf, vjust = -0.8,
            size = 2.7, colour = "grey40") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(title = "Disparity and optimality through the Mesozoic",
       subtitle = paste("Pareto rank ratio replaces the former averaged score.",
                        CAPTION_CI),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())

ggsave("output/plots/TEMPORAL_03_disparity_optimality.png",
       p_disp, width = 11, height = 4, dpi = 300)
ggsave("output/plots_PDF/TEMPORAL_03_disparity_optimality.pdf",
       p_disp, width = 11, height = 4)
cat("Written: TEMPORAL_03_disparity_optimality\n\n")
