###################################################################################
# 2. Body size
###################################################################################
# Size as its own subject rather than an unacknowledged driver of the shape
# metrics. If size rises across the Cretaceous while shape traits stay flat,
# that separation is itself the result.
size_vars <- intersect(c("Wingspan_cm", "Mass_kg", "Wing_area_cm2"), present)

p_size <- make_time_panel(size_vars, "Body size through the Mesozoic",
                          paste("Log10 scale.", CAPTION))
save_fig(p_size, "TEMPORAL_02_body_size", length(size_vars))

###################################################################################
# 3. Disparity and optimality
###################################################################################
disp_long <- rbind(
  data.frame(Time_Bin = disparity_metrics$Time_Bin, n = disparity_metrics$n,
             variable = "sov", label = "Shape disparity (sum of variance)",
             mean = disparity_metrics$sov,
             lo = disparity_metrics$sov_lo, hi = disparity_metrics$sov_hi),
  data.frame(Time_Bin = disparity_metrics$Time_Bin, n = disparity_metrics$n,
             variable = "mpd", label = "Morphospace expansion (mean pairwise distance)",
             mean = disparity_metrics$mpd,
             lo = disparity_metrics$mpd_lo, hi = disparity_metrics$mpd_hi))

if ("pareto_rank_ratio" %in% present) {
  prr <- temporal_metrics[temporal_metrics$variable == "pareto_rank_ratio", ]
  disp_long <- rbind(disp_long,
                     data.frame(Time_Bin = prr$Time_Bin, n = prr$n,
                                variable = "prr", label = "Pareto rank ratio",
                                mean = prr$mean, lo = prr$lo, hi = prr$hi))
}

DISP_COLS <- c("Shape disparity (sum of variance)" = "#1F6FB4",
               "Morphospace expansion (mean pairwise distance)" = "#B5651D",
               "Pareto rank ratio" = "#1b7837")
disp_long$label <- factor(disp_long$label, levels = names(DISP_COLS))

p_disp <- ggplot(disp_long, aes(x = Time_Bin, y = mean, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = label), alpha = 0.18, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi, colour = label),
                width = 0.16, linewidth = 0.55) +
  geom_line(aes(colour = label), linewidth = 0.85) +
  geom_point(aes(colour = label), size = 2.7) +
  geom_text(aes(label = paste0("n=", n)), y = -Inf, vjust = -0.7,
            size = 2.5, colour = "grey35") +
  scale_colour_manual(values = DISP_COLS, guide = "none") +
  scale_fill_manual(values = DISP_COLS, guide = "none") +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(title = "Disparity and optimality through the Mesozoic",
       subtitle = paste("Pareto rank ratio replaces the former averaged score.",
                        CAPTION), x = NULL, y = NULL) +
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

save_fig(p_disp, "TEMPORAL_03_disparity_optimality", nlevels(disp_long$label))

cat("\nBefore interpreting any trend, check the PGLS in Msc_phylo_signal_v3.R:\n")
