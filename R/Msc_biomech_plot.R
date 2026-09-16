cat("\n================================================================================\n")
cat("BIOMECHANICAL MORPHOSPACE PANELS\n")
cat("================================================================================\n\n")
# stress_index removed - von_mises_stress is the only stress metric.
METRIC_LABELS <- c(
  von_mises_stress  = "Von Mises stress\n(scale-invariant)",
  r2_hat            = "Second moment\nof area",
  aspect_ratio      = "Aspect ratio",
  wing_loading_ratio      = "Wing loading\n(N/m\u00b2)",
  wing_curvature    = "Wing curvature",
  shape_complexity  = "Shape complexity",
  pareto_rank_ratio = "Pareto\nrank ratio")

surf_label <- if (BACKGROUND_VAR %in% names(METRIC_LABELS)) {
  METRIC_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR

cat("Surface variable:", BACKGROUND_VAR, "| colour bar reads:", gsub("\n", " ", surf_label), "\n")

# --------------------------------------------------------------------------------
# 2. Convex hulls, rebuilt from the points actually plotted and closed
# --------------------------------------------------------------------------------
build_hulls <- function(group_col = "clade") {
  if (!group_col %in% colnames(performance_data_clean)) return(NULL)
  d <- performance_data_clean[!is.na(performance_data_clean[[group_col]]) & is.finite(performance_data_clean$PC1) &is.finite(performance_data_clean$PC2), ]
  if (nrow(d) < 3) return(NULL)
  do.call(rbind, lapply(unique(d[[group_col]]), function(g) {
    dg <- d[d[[group_col]] == g, ]
    if (nrow(dg) < 3) return(NULL)
    h <- chull(dg$PC1, dg$PC2); h <- c(h, h[1])
    data.frame(PC1 = dg$PC1[h], PC2 = dg$PC2[h], grp = g)
  }))
}

hulls_clade <- build_hulls("clade")
has_hulls <- !is.null(hulls_clade) && nrow(hulls_clade) > 0
cat(sprintf("Clade hulls: %d groups\n\n",
            if (has_hulls) length(unique(hulls_clade$grp)) else 0))

ALL_CLADES <- if ("clade" %in% colnames(performance_data_clean)) {
  unique(as.character(performance_data_clean$clade))
} else character(0)
CLADE_COLS_MATCHED <- match_clade_colour(ALL_CLADES)
unmapped_cl <- ALL_CLADES[is.na(CLADE_COLS_MATCHED)]
if (length(unmapped_cl) > 0) {
  cat("NOTE: no CLADE_COLS entry for:", paste(unmapped_cl, collapse = ", "),
      "- these get a grey border.\n")
  CLADE_COLS_MATCHED[is.na(CLADE_COLS_MATCHED)] <- "grey50"
}
cat("Clade colour mapping in use:\n")
print(CLADE_COLS_MATCHED)
cat("\n")

# --------------------------------------------------------------------------------
# 3. Panel builder
# --------------------------------------------------------------------------------
make_biomech_panel <- function(colour_var, panel_title, palette = "plasma",
                               log_colour = FALSE, use_stress_palette = FALSE,
                               hollow_points = FALSE) 
{
  if (!colour_var %in% colnames(performance_data_clean)) {
    cat("SKIPPED", panel_title, "- column", colour_var, "does not exist.\n")
    cat("  Present:", paste(intersect(names(METRIC_LABELS), colnames(performance_data_clean)), collapse = ", "), "\n")
    return(NULL)
  }
  
  d <- performance_data_clean
  v <- d[[colour_var]]
  
  d$.clade <- if ("clade" %in% colnames(d)) as.character(d$clade) else NA_character_
  
  if (log_colour && any(v <= 0, na.rm = TRUE)) {
    cat(sprintf("  %s: %d non-positive values - linear scale used.\n",colour_var, sum(v <= 0, na.rm = TRUE)))
    log_colour <- FALSE
  }
  d$.col <- if (log_colour) log(v) else v
  d <- d[is.finite(d$.col), ]
  
  if (nrow(d) < 5) {
    cat("SKIPPED", panel_title, "- only", nrow(d), "finite values.\n"); return(NULL)
  }
  
  pt_label_raw <- if (colour_var %in% names(METRIC_LABELS)) {
    METRIC_LABELS[[colour_var]]
  } else colour_var
  if (log_colour) pt_label_raw <- paste0("log(", gsub("\n", " ", pt_label_raw), ")")
  pt_label <- paste0("Point fill - ", gsub("\n", " ", pt_label_raw))
  bg_label <- paste0("Background - ", gsub("\n", " ", surf_label))
  
  cat(sprintf("  %s: %d specimens | legend: %s | points: %s\n",
              colour_var, nrow(d), gsub("\n", " ", pt_label_raw),
              if (hollow_points) "hollow (clade border only)" else "filled"))
  
  p <- ggplot() +geom_tile(data = stress_grid, aes(x = PC1, y = PC2, fill = stress),alpha = 0.55) +
    scale_fill_stress(name = bg_label, palette = "dark")
  
  if (!is.null(wing_bg)) {
    p <- p + geom_polygon(data = wing_bg, aes(x = PC1, y = PC2, group = grid_id), fill = "#2a2a2a", colour = "#1a1a1a",alpha = 0.18, linewidth = 0.25)
  }
  
  p <- p + ggnewscale::new_scale_colour()
  
  if (has_hulls) {
    p <- p + geom_path(data = hulls_clade, aes(x = PC1, y = PC2, group = grp, colour = grp),
                       linewidth = 1.0, alpha = 0.9, show.legend = FALSE)
  }
  
  if (hollow_points) {
    p <- p +
      geom_point(data = d, aes(x = PC1, y = PC2, colour = .clade),
                 shape = 21, fill = NA, size = 3.2, stroke = 1.0, alpha = 0.9) +
      scale_colour_manual(values = CLADE_COLS_MATCHED, na.value = "grey50",
                          name = "Point border - Clade")
  } else {
    p <- p + ggnewscale::new_scale_fill()
    
    fill_scale <- if (use_stress_palette) {
      ggplot2::scale_fill_gradientn(colors = STRESS_PALETTES$dark, name = pt_label, na.value = "white")
    } else {
      scale_fill_viridis_c(option = palette, name = pt_label)
    }
    
    p <- p +
      geom_point(data = d, aes(x = PC1, y = PC2, fill = .col, colour = .clade),
                 shape = 21, size = 3.2, stroke = 0.9, alpha = 0.9) +
      fill_scale +
      scale_colour_manual(values = CLADE_COLS_MATCHED, na.value = "grey50",
                          name = "Point border - Clade")
  }
  
  p +
    labs(title = panel_title,
         subtitle = paste0("Surface: ", gsub("\n", " ", surf_label),
                           ". Dark outlines: wing profiles across the morphospace.",
                           if (hollow_points) " Points: clade border only, no fill." else ""),
         x = pc_lab(1), y = pc_lab(2)) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0, size = 9, colour = "grey35"),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey90"),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6))
}

# --------------------------------------------------------------------------------
# 4. The three panels
# --------------------------------------------------------------------------------
# No fallback: von_mises_stress must exist, or this stops with a clear message
# instead of silently substituting a different metric.
if (!"von_mises_stress" %in% colnames(performance_data_clean)) {
  stop("von_mises_stress not found in performance_data_clean. ",
       "Source Msc_performance_metrics_v3.R first.")
}
stress_col <- "von_mises_stress"

PANELS <- list(
  list(var = "aspect_ratio", title = "(a) Aspect ratio", pal = "plasma", log = TRUE, stress = FALSE, hollow = FALSE, file = "CLADE_01A_aspect_ratio"),
  list(var = "r2_hat", title = "(b) Second moment of area", pal = "viridis", log = FALSE, stress = FALSE, hollow = FALSE, file = "CLADE_01B_second_moment"),
  list(var = stress_col, title = "(c) Von Mises stress", pal = "mako", log = TRUE, stress = TRUE, hollow = TRUE, file = "CLADE_01C_von_mises_stress"))

plots <- list()
for (pn in PANELS) {
  if (is.na(pn$var)) next
  p <- make_biomech_panel(pn$var, pn$title, pn$pal, pn$log, pn$stress, pn$hollow)
  if (is.null(p)) next
  plots[[pn$file]] <- p
  ggsave(sprintf("output/plots_PDF/biomechanic_plot/%s.pdf", pn$file),
         p, width = 18, height = 11)
  cat("Written:", pn$file, "\n")
}