# Make sure the plot folders and input table are ready.
for (dd in c("output/plots", "output/plots_PDF")) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("performance_data_clean")) stop("performance_data_clean not found.")
if (!exists("PTERODACT")) {
  PTERODACT <- c("Azhdarchoidea", "Ctenochasmatoidea", "Dsungaripteroidea",
                 "Ornithocheiromorpha", "Pteranodontia", "basal pterodactyloidea")
}

# Same convention as Msc_phylo_v3_1.R's order_group: Darwinoptera on the
# Non-Pterodactyliform side.
perf_pareto_df <- performance_data_clean
perf_pareto_df$panel_group <- ifelse(
  perf_pareto_df$clade %in% PTERODACT & perf_pareto_df$clade != "Darwinoptera",
  "Pterodactyliform", "Non-Pterodactyliform")

# The three genuine Pareto objectives - not second_moment_DIAGNOSTIC, which
# enters the von_mises_stress formula as a divisor and would redraw an
# algebraic identity rather than an independent comparison.
METRICS_D <- c("aspect_ratio", "r2_hat", "von_mises_stress")
missing_d <- setdiff(METRICS_D, names(perf_pareto_df))
if (length(missing_d) > 0) {
  stop("Missing column(s): ", paste(missing_d, collapse = ", "),
       "\n  Run Msc_performance_metrics_v3.R first.")
}

d <- perf_pareto_df[is.finite(perf_pareto_df$aspect_ratio) & perf_pareto_df$aspect_ratio > 0 &
                      is.finite(perf_pareto_df$r2_hat) & perf_pareto_df$r2_hat > 0 &
                      is.finite(perf_pareto_df$von_mises_stress) & perf_pareto_df$von_mises_stress > 0 &
                      !is.na(perf_pareto_df$panel_group), ]
cat(sprintf("Specimens with all three metrics finite and positive: %d / %d\n",
            nrow(d), nrow(perf_pareto_df)))
print(table(d$panel_group))
cat("\n")

GROUP_COLS <- c("Non-Pterodactyliform" = "#A64CA6", "Pterodactyliform" = "#F2C230")

has_theoretical <- exists("theoretical_data") &&
  all(c("self_intersecting") %in% names(theoretical_data))

if (!has_theoretical) {
  cat("theoretical_data not found - envelope will show real specimens' own hull only.\n")
  cat("  Source Msc_theoretical_spaces_v3.R first to include the theoretical shapes too.\n\n")
}

# --------------------------------------------------------------------------------
# 2D Pareto front (skyline), computed on the POOLED data - group membership
# plays no part in deciding who is on the front, only in colouring points
# afterwards. Returns a logical vector aligned with the INPUT order.
# --------------------------------------------------------------------------------
pareto_front_2d <- function(x, y, x_max = TRUE, y_max = TRUE) {
  n <- length(x)
  keep <- rep(FALSE, n)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) return(keep)

  idx <- which(ok)
  xo <- x[idx]; yo <- y[idx]
  ord <- order(xo, decreasing = x_max)
  xo <- xo[ord]; yo <- yo[ord]; idx <- idx[ord]

  best <- if (y_max) -Inf else Inf
  for (i in seq_along(xo)) {
    better <- if (y_max) yo[i] > best else yo[i] < best
    if (better) { keep[idx[i]] <- TRUE; best <- yo[i] }
  }
  keep
}

##################################################################################
# Build the outline of the possible solution space.
##################################################################################
build_theoretical_envelope <- function(xcol, ycol) {
  pts <- data.frame(x = numeric(0), y = numeric(0))

  if (has_theoretical && all(c(xcol, ycol) %in% names(theoretical_data))) {
    td <- theoretical_data[!theoretical_data$self_intersecting &
                             is.finite(theoretical_data[[xcol]]) & theoretical_data[[xcol]] > 0 &
                             is.finite(theoretical_data[[ycol]]) & theoretical_data[[ycol]] > 0, ]
    if (nrow(td) > 0) pts <- rbind(pts, data.frame(x = td[[xcol]], y = td[[ycol]]))
  }

  real_ok <- is.finite(d[[xcol]]) & d[[xcol]] > 0 & is.finite(d[[ycol]]) & d[[ycol]] > 0
  pts <- rbind(pts, data.frame(x = d[[xcol]][real_ok], y = d[[ycol]][real_ok]))

  if (nrow(pts) < 3) return(NULL)
  h <- chull(pts$x, pts$y); h <- c(h, h[1])
  data.frame(x = pts$x[h], y = pts$y[h])
}

################################################################################
# Build one Pareto panel for a chosen pair of measures.
################################################################################
make_perf_panel <- function(xcol, ycol, xlab, ylab, x_max, y_max, letter,
                            log_x = FALSE, log_y = FALSE) {
  on_front <- pareto_front_2d(d[[xcol]], d[[ycol]], x_max, y_max)
  front_pts <- d[on_front, c(xcol, ycol)]
  front_pts <- front_pts[order(front_pts[[xcol]]), ]

  front_tab <- table(d$panel_group[on_front])
  cat(sprintf("  %s vs %s: %d specimens on the common front (%s)\n",
              xcol, ycol, nrow(front_pts),
              paste(sprintf("%s: %d", names(front_tab), front_tab), collapse = ", ")))

  envelope <- build_theoretical_envelope(xcol, ycol)

  p <- ggplot()

  if (!is.null(envelope)) {
    p <- p +
      geom_polygon(data = envelope, aes(x = x, y = y, fill = "Theoretical aerodynamic solutions"),
                   colour = NA, alpha = 0.35) +
      scale_fill_manual(name = NULL, values = c("Theoretical aerodynamic solutions" = "grey55")) +
      ggnewscale::new_scale_fill()
  }

  p <- p +
    geom_point(data = d, aes(x = .data[[xcol]], y = .data[[ycol]], fill = panel_group),
               shape = 21, colour = "black", stroke = 0.3, size = 2.3, alpha = 0.9) +
    scale_fill_manual(name = NULL, values = GROUP_COLS)
  p <- p + if (log_x) scale_x_log10(labels = scales::label_log()) else scale_x_continuous()
  p <- p + if (log_y) scale_y_log10(labels = scales::label_log()) else scale_y_continuous()
  if (nrow(front_pts) >= 2) {
    p <- p + geom_line(data = front_pts, aes(x = .data[[xcol]], y = .data[[ycol]], colour = "Pareto front"),
                       linewidth = 0.9, inherit.aes = FALSE) +
      scale_colour_manual(name = NULL, values = c("Pareto front" = "black"))
  }

  p +
    labs(subtitle = sprintf("(%s)", letter), x = xlab, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "grey20", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          legend.position = "none",
          plot.subtitle = element_text(face = "bold", size = 11))
}

cat("Common Pareto front composition per panel:\n")
p1 <- make_perf_panel("aspect_ratio", "r2_hat","Aspect Ratio", "Second moment of area (r2_hat)",x_max = TRUE, y_max = TRUE, letter = "a", log_x = FALSE, log_y = FALSE)

p2 <- make_perf_panel("aspect_ratio", "von_mises_stress",
                      "Aspect Ratio", "Von Mises stress (scale-invariant, log)",
                      x_max = TRUE, y_max = FALSE, letter = "b",
                      log_x = FALSE, log_y = TRUE)

p3 <- make_perf_panel("von_mises_stress", "r2_hat",
                      "Von Mises stress (scale-invariant, log)", "Second moment of area (r2_hat)",
                      x_max = FALSE, y_max = TRUE, letter = "c",
                      log_x = TRUE, log_y = FALSE)
cat("\n")

# Build one shared legend from panel b so all guides stay together.
legend_src <- p2 + theme(legend.position = "bottom", legend.box = "horizontal")
shared_legend <- cowplot::get_legend(legend_src)

combined_d <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
combined_d <- cowplot::plot_grid(combined_d, shared_legend, ncol = 1, rel_heights = c(1, 0.12))

title_d <- cowplot::ggdraw() +cowplot::draw_label( "Pairwise optimality comparison", fontface = "bold", size = 15, x = 0, hjust = 0)

final_d <- cowplot::plot_grid(title_d, combined_d, ncol = 1, rel_heights = c(0.1, 1))
ggsave("output/plots_PDF/CLADE_15D_pairwise_performance_pareto.pdf",final_d, width = 15, height = 6.5)
cat("Written: CLADE_15D_pairwise_performance_pareto\n\n")
