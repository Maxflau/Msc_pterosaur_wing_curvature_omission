# --------------------------------------------------------------------------------
# TAPH_01 / TAPH_02 — draw the preservation checks used in the supplement.
# --------------------------------------------------------------------------------
if (exists("POSTCRANIAL") && "n_elements" %in% colnames(POSTCRANIAL)) {
  tp <- POSTCRANIAL
  tp$completeness <- tp$n_elements / length(GM_ELEMENTS)
  tp[[ENV_VAR]] <- trimws(as.character(tp[[ENV_VAR]]))
  tp <- tp[!is.na(tp[[ENV_VAR]]) & tp[[ENV_VAR]] != "" & is.finite(tp$completeness), ]
  tp$setting <- tp[[ENV_VAR]]
  tp <- tp[tp$setting %in% names(which(table(tp$setting) >= MIN_N)), ]

  p_t1 <- ggplot(tp, aes(x = reorder(setting, completeness, FUN = median),
                         y = completeness, fill = setting)) +
    geom_violin(alpha = 0.45, scale = "width", colour = "grey35", linewidth = 0.3) +
    geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.08, size = 0.8, alpha = 0.5, colour = "grey25") +
    scale_fill_viridis_d(option = "rocket", guide = "none", begin = 0.25, end = 0.85) +
    coord_flip() +
    labs(title = "Skeletal completeness by depositional setting",
         subtitle = "A difference here means settings preserve differently, so any metric difference between them is confounded with preservation.",
         x = NULL, y = "Fraction of postcranial elements measured") +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey35"),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))

  ggsave("output/plots/TAPH_01_completeness.png", p_t1,
         width = 10, height = 6, dpi = 300)
  ggsave("output/plots_PDF/TAPH_01_completeness.pdf", p_t1,
         width = 10, height = 6)
  cat("Written: TAPH_01_completeness\n")
  # Turn the bias check into a long table so each metric can get its own panel.
  bias_long <- do.call(rbind, lapply(intersect(names(LABELS), colnames(tp)),
                                     function(m) {
                                       data.frame(completeness = tp$completeness, metric = LABELS[[m]],
                                                  value = tp[[m]], stringsAsFactors = FALSE)
                                     }))
  bias_long <- bias_long[is.finite(bias_long$value) &
                           is.finite(bias_long$completeness), ]

  p_t2 <- ggplot(bias_long, aes(x = completeness, y = value)) +
    geom_point(size = 0.9, alpha = 0.45, colour = "grey30") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#b2182b",
                fill = "#b2182b", alpha = 0.15, linewidth = 0.7) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    labs(title = "Preservation bias: metrics against skeletal completeness",
         subtitle = "A visible slope means part of that metric's variance is preservation, not anatomy. See S20 for the correlations.",
         x = "Fraction of postcranial elements measured", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey95", colour = "grey40"),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
          plot.subtitle = element_text(size = 8, colour = "grey35"))

  ggsave("output/plots/TAPH_02_completeness_vs_metrics.png", p_t2,
         width = 12, height = 7, dpi = 300)
  ggsave("output/plots_PDF/TAPH_02_completeness_vs_metrics.pdf", p_t2,
         width = 12, height = 7)
  cat("Written: TAPH_02_completeness_vs_metrics\n")
} else {
  cat("POSTCRANIAL not found - taphonomy figures skipped.\n")
  cat("Source Msc_stats_06a_postcranial.R first.\n")
}
cat("\n")
