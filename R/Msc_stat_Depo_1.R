cat("\n--- ENVIRONMENT AND TAPHONOMY FIGURES ---\n\n")

library(ggplot2)
for (dd in c("output/plots", "output/plots_PDF")) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("METRICS")) stop("Source Msc_stats_00_setup.R first.")

LABELS <- c(aspect_ratio = "Aspect ratio", r2_hat = "Second moment of area",
            stress_root = "Root bending stress index",
            wing_curvature = "Wing curvature", shape_complexity = "Shape complexity",
            wing_loading = "Wing loading (N/m2)",
            pareto_rank_ratio = "Pareto rank ratio")

ENV_VAR <- intersect(c("Depositional.settings.paleoenvironment", "Environment"),
                     colnames(performance_data_clean))[1]
if (is.na(ENV_VAR)) stop("No environmental column found.")
cat("Environmental factor:", ENV_VAR, "\n")

d <- performance_data_clean
d[[ENV_VAR]] <- trimws(as.character(d[[ENV_VAR]]))
d <- d[!is.na(d[[ENV_VAR]]) & d[[ENV_VAR]] != "", ]
keep_lv <- names(which(table(d[[ENV_VAR]]) >= MIN_N))
d <- d[d[[ENV_VAR]] %in% keep_lv, ]
d$setting <- d[[ENV_VAR]]
cat(sprintf("Settings with n >= %d: %d | specimens: %d\n\n", MIN_N,
            length(keep_lv), nrow(d)))

# --------------------------------------------------------------------------------
# ENV_01 — metrics by setting
# --------------------------------------------------------------------------------
PLOT_METRICS <- intersect(names(LABELS), colnames(d))

long <- do.call(rbind, lapply(PLOT_METRICS, function(m) {
  data.frame(setting = d$setting, metric = LABELS[[m]], value = d[[m]],
             stringsAsFactors = FALSE)
}))
long <- long[is.finite(long$value), ]
long$metric <- factor(long$metric, levels = LABELS[PLOT_METRICS])

n_lab <- as.data.frame(table(d$setting))
names(n_lab) <- c("setting", "n")
n_lab$metric <- factor(LABELS[PLOT_METRICS[1]], levels = levels(long$metric))

p_env <- ggplot(long, aes(x = setting, y = value, fill = setting)) +
  geom_violin(alpha = 0.45, scale = "width", colour = "grey35", linewidth = 0.3) +
  geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.08, size = 0.6, alpha = 0.4, colour = "grey25") +
  geom_text(data = n_lab, aes(x = setting, label = paste0("n=", n)), y = -Inf,
            vjust = -0.6, size = 2.3, colour = "grey40", inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_viridis_d(option = "cividis", guide = "none", begin = 0.2, end = 0.9) +
  labs(title = paste("Performance metrics by", ENV_VAR),
       subtitle = "Read with S14: settings are not independent of geological period.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_rect(fill = "grey95", colour = "grey40"),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
        plot.title = element_text(face = "bold", size = 13))

h <- 3.3 * ceiling(length(PLOT_METRICS) / 2)
ggsave("output/plots/ENV_01_metrics_by_setting.png", p_env,
       width = 11, height = h, dpi = 300, limitsize = FALSE)
ggsave("output/plots_PDF/ENV_01_metrics_by_setting.pdf", p_env,
       width = 11, height = h, limitsize = FALSE)
cat("Written: ENV_01_metrics_by_setting\n")

# --------------------------------------------------------------------------------
# ENV_02 — morphospace occupation by setting
# --------------------------------------------------------------------------------
if (all(c("PC1", "PC2") %in% colnames(d))) {
  
  hulls <- do.call(rbind, lapply(unique(d$setting), function(s) {
    ds <- d[d$setting == s & is.finite(d$PC1) & is.finite(d$PC2), ]
    if (nrow(ds) < 3) return(NULL)
    h <- chull(ds$PC1, ds$PC2); h <- c(h, h[1])   # close the ring
    data.frame(PC1 = ds$PC1[h], PC2 = ds$PC2[h], setting = s)
  }))
  
  p_ms <- ggplot() +
    geom_point(data = performance_data_clean, aes(x = PC1, y = PC2),
               colour = "grey88", size = 0.8) +
    {if (!is.null(hulls))
      geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = setting,
                                     fill = setting, colour = setting),
                   alpha = 0.25, linewidth = 0.6)} +
    geom_point(data = d, aes(x = PC1, y = PC2, colour = setting), size = 2.2) +
    scale_fill_viridis_d(option = "cividis", name = NULL, begin = 0.2, end = 0.9) +
    scale_colour_viridis_d(option = "cividis", name = NULL, begin = 0.2, end = 0.9) +
    labs(title = paste("Performance morphospace by", ENV_VAR),
         subtitle = "Grey: full sample.",
         x = sprintf("PC1 (%.1f%%)", summary(pca_performance)$importance[2, 1] * 100),
         y = sprintf("PC2 (%.1f%%)", summary(pca_performance)$importance[2, 2] * 100)) +
    theme_bw(base_size = 12) +
    theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 13))
  
  ggsave("output/plots/ENV_02_morphospace_by_setting.png", p_ms,
         width = 10, height = 7, dpi = 300)
  ggsave("output/plots_PDF/ENV_02_morphospace_by_setting.pdf", p_ms,
         width = 10, height = 7)
  cat("Written: ENV_02_morphospace_by_setting\n")
}