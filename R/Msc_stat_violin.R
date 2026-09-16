if (!exists("METRICS")) stop("Source Msc_stats_00_setup.R first.")

for (d in c("output/plots", "output/plots_PDF")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

VIOLIN_LEVELS <- intersect(c("clade", "Order", "Diet_combined", "Diet.1"),
                           colnames(performance_data_clean))
VIOLIN_METRICS <- intersect(c("aspect_ratio", "r2_hat", "stress_root",
                              "wing_curvature", "shape_complexity",
                              "wing_loading", "pareto_rank_ratio"),
                            colnames(performance_data_clean))

METRIC_LABELS <- c(aspect_ratio = "Aspect ratio", r2_hat = "Second moment of area",
                   stress_root = "Root bending stress index",
                   wing_curvature = "Wing curvature", shape_complexity = "Shape complexity",
                   wing_loading = "Wing loading (N/m2)",
                   pareto_rank_ratio = "Pareto rank ratio")

cat("Levels:", paste(VIOLIN_LEVELS, collapse = ", "), "\n")
cat("Metrics:", paste(VIOLIN_METRICS, collapse = ", "), "\n\n")

# --------------------------------------------------------------------------------
# Builder: one faceted figure per grouping level
# --------------------------------------------------------------------------------
make_violin <- function(group_col) {
  
  d <- performance_data_clean[!is.na(performance_data_clean[[group_col]]) &
                                performance_data_clean[[group_col]] != "", ]
  d$.grp <- as.character(d[[group_col]])
  
  # Groups of one or two are kept as points but excluded from the violin layer
  n_by <- table(d$.grp)
  big  <- names(n_by)[n_by >= MIN_N]
  if (length(big) < 2) { cat("SKIPPED", group_col, "- too few usable groups\n"); return(NULL) }
  
  long <- do.call(rbind, lapply(VIOLIN_METRICS, function(m) {
    data.frame(grp = d$.grp, metric = METRIC_LABELS[m], value = d[[m]],
               stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$value), ]
  long$metric <- factor(long$metric, levels = METRIC_LABELS[VIOLIN_METRICS])
  
  # Order groups by median of the first metric, so the panels read consistently
  ref <- long[long$metric == METRIC_LABELS[VIOLIN_METRICS[1]], ]
  ord <- names(sort(tapply(ref$value, ref$grp, median, na.rm = TRUE)))
  long$grp <- factor(long$grp, levels = ord)
  
  n_lab <- data.frame(grp = factor(names(n_by), levels = ord),
                      n = as.integer(n_by),
                      metric = factor(METRIC_LABELS[VIOLIN_METRICS[1]],
                                      levels = levels(long$metric)))
  
  cat(sprintf("  %s: %d groups (%d with n >= %d)\n",
              group_col, length(n_by), length(big), MIN_N))
  
  ggplot(long, aes(x = grp, y = value)) +
    geom_violin(data = long[long$grp %in% big, ],
                aes(fill = grp), alpha = 0.45, colour = "grey35",
                linewidth = 0.35, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.85,
                 outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.09, size = 0.7, alpha = 0.45, colour = "grey25") +
    geom_text(data = n_lab, aes(x = grp, label = paste0("n=", n)),
              y = -Inf, vjust = -0.6, size = 2.3, colour = "grey40",
              inherit.aes = FALSE) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_fill_viridis_d(option = "mako", guide = "none", begin = 0.15, end = 0.9) +
    labs(title = paste("Performance metrics by", group_col),
         subtitle = paste("Violins only for groups with n >=", MIN_N,
                          "- smaller groups show points alone."),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey95", colour = "grey40"),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 13))
}

# --------------------------------------------------------------------------------
# Build and write
# --------------------------------------------------------------------------------
for (g in VIOLIN_LEVELS) {
  p <- make_violin(g)
  if (is.null(p)) next
  n_g <- length(unique(performance_data_clean[[g]]))
  w <- max(10, min(20, 5 + n_g * 0.45))
  h <- 3.4 * ceiling(length(VIOLIN_METRICS) / 2)
  nm <- paste0("VIOLIN_", tolower(gsub("\\.", "_", g)))
  ggsave(sprintf("output/plots/%s.png", nm), p,
         width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(sprintf("output/plots_PDF/%s.pdf", nm), p,
         width = w, height = h, limitsize = FALSE)
  cat(sprintf("Written: %s (%.1f x %.1f in)\n", nm, w, h))
}