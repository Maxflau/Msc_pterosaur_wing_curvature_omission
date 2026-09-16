# ── 4. BASE THEME ─────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title        = element_text(face="bold", size=12),
    plot.subtitle     = element_text(size=9, colour="grey40"),
    plot.margin       = unit(c(4, 105, 4, 4), "pt"),
    axis.title        = element_text(size=10),
    legend.position   = "right",
    legend.title      = element_text(size=8.5, face="bold"),
    legend.text       = element_text(size=7.5),
    legend.key.height = unit(0.42, "cm"),
    legend.key.width  = unit(0.42, "cm"),
    legend.spacing.y  = unit(0.12, "cm"),
    legend.box        = "vertical",
    panel.border      = element_rect(colour="grey30", fill=NA, linewidth=0.5),
    panel.grid.major  = element_line(colour="grey88", linewidth=0.5),
    panel.grid.minor  = element_blank())

# Metrics spanning more than an order of magnitude; everything else is a ratio
LOG_METRICS <- c("von_mises_stress", "wing_loading_ratio")

# ── 5. CORE PLOT FUNCTION ─────────────────────────────────────────────────────
# border (colour) = group | interior (fill) = metric value
make_metric_plot <- function(
    metric_col, metric_label,
    metric_log      = FALSE,
    metric_pal      = "viridis",
    metric_rev      = FALSE,
    metric_breaks   = NULL,
    group_var       = "clade",
    group_cols      = CLADE_COLS,
    hull_fill       = FALSE,
    hull_alpha      = 0.15,
    clade_colour    = FALSE,
    transparent_pts = FALSE,
    dark_bg         = FALSE,
    legend_bottom   = FALSE,
    perm_factor     = "clade",
    title           = NULL,
    subtitle        = NULL) {
  
  if (!metric_col %in% names(shape_perf)) {
    message("Column not found: ", metric_col); return(NULL)
  }
  
  mvals <- as.numeric(shape_perf[[metric_col]])
  # Non-positive values are dropped, not set to 1e-9: a fabricated near-zero
  # becomes an extreme outlier once logged and stretches the whole colour ramp
  if (metric_log && metric_col %in% LOG_METRICS) {
    mvals[!is.na(mvals) & mvals <= 0] <- NA
    mvals <- log10(mvals)
  } else if (metric_log) {
    message(metric_col, " is a ratio - log ignored.")
  }
  
  pt_df <- shape_perf %>% mutate(mval = mvals) %>% filter(is.finite(mval))
  if (nrow(pt_df) < 5) { message("Too few finite values: ", metric_col); return(NULL) }
  
  hulls2 <- make_hulls(shape_perf, group_var)
  pv <- get_perm(perm_factor)
  
  VIRIDIS_OPTIONS <- c("magma","inferno","plasma","viridis","cividis",
                       "rocket","mako","turbo")
  use_viridis <- metric_pal %in% VIRIDIS_OPTIONS
  
  # Breaks outside the data range are silently dropped, so a hard-coded list
  # that no longer matches will desynchronise from its labels
  if (!is.null(metric_breaks)) {
    r <- range(pt_df$mval, na.rm = TRUE)
    metric_breaks <- metric_breaks[metric_breaks >= r[1] & metric_breaks <= r[2]]
    if (length(metric_breaks) < 2) metric_breaks <- NULL
  }
  
  metric_scale <- if (use_viridis) {
    scale_fill_viridis_c(option = metric_pal,
                         direction = if (metric_rev) -1 else 1,
                         name = metric_label, breaks = metric_breaks,
                         guide = guide_colorbar(order=2, barheight=unit(3.5,"cm"),
                                                barwidth=unit(0.32,"cm"),
                                                title.theme=element_text(size=7, face="bold")))
  } else {
    scale_fill_distiller(palette = metric_pal,
                         direction = if (metric_rev) -1 else 1,
                         name = metric_label, breaks = metric_breaks,
                         guide = guide_colorbar(order=2, barheight=unit(3.5,"cm"),
                                                barwidth=unit(0.32,"cm"),
                                                title.theme=element_text(size=7, face="bold")))
  }
  
  p <- ggplot()
  
  # 1 ── Stress background: breaks left to ggplot, only the format is fixed
  if (!is.null(stress_bg)) {
    p <- p + geom_raster(data=stress_bg, aes(x, y, fill=z), interpolate=TRUE)
    p <- p + if (dark_bg) {
      scale_fill_gradientn(
        colours = c("#FFFFFF","#DEEBF7","#9ECAE1","#4292C6","#08519C","#08306B"),
        name = "Von Mises stress\n(log10)",
        labels = scales::label_number(accuracy = 0.01),
        guide = guide_colorbar(order=99, barheight=unit(4,"cm"),
                               barwidth=unit(0.32,"cm"),
                               title.theme=element_text(size=7, face="bold")))
    } else {
      scale_fill_distiller(palette="Blues", direction=1,
                           name = "Von Mises stress\n(log10)",
                           labels = scales::label_number(accuracy = 0.01),
                           guide = guide_colorbar(order=99, barheight=unit(3.5,"cm"),
                                                  barwidth=unit(0.32,"cm"),
                                                  title.theme=element_text(size=7, face="bold")))
    }
  }
  
  # 2 ── Region outside the occupied envelope, if the pipeline supplied one
  if (nrow(impossible_region) > 0) {
    p <- p + geom_tile(data=impossible_region, aes(x, y),
                       fill="grey80", alpha=0.75, inherit.aes=FALSE)
  }
  
  # 3 ── EFA grid wings
  p <- p + geom_polygon(data=grid_dense, aes(gx, gy, group=grid_id),
                        fill="white", colour="grey40", linewidth=0.07, alpha=0.45)
  
  # 4 + 5 ── Hulls and points; new_scale_fill() separates background from metric
  if (clade_colour) {
    p <- p +
      geom_polygon(data=hulls2, aes(shapePC1, shapePC2, group=grp, colour=grp),
                   fill=NA, linewidth=0.72) +
      scale_colour_manual(values=group_cols, name="Clade",
                          guide=guide_legend(order=1,
                                             override.aes=list(shape=21, fill="grey70",
                                                               size=3.5, stroke=1.0, linewidth=NA)))
    if (transparent_pts) {
      p <- p + new_scale_fill() +
        geom_point(data=pt_df, aes(shapePC1, shapePC2, colour=.data[[group_var]]),
                   shape=21, fill=NA, size=3.0, stroke=0.75)
    } else {
      p <- p + new_scale_fill() +
        geom_point(data=pt_df, aes(shapePC1, shapePC2, fill=mval,
                                   colour=.data[[group_var]]),
                   shape=21, size=3.0, stroke=0.60) +
        metric_scale +
        guides(fill = guide_colorbar(order=2, barheight=unit(4,"cm"),
                                     barwidth=unit(0.38,"cm"),
                                     title.theme=element_text(size=8, face="bold")))
    }
  } else {
    p <- p +
      geom_polygon(data=hulls2, aes(shapePC1, shapePC2, group=grp),
                   fill=NA, colour="black", linewidth=0.65) +
      new_scale_fill() +
      geom_point(data=pt_df, aes(shapePC1, shapePC2, fill=mval),
                 shape=21, colour="grey15", size=2.6, stroke=0.45) +
      metric_scale
  }
  
  plot_title <- if (!is.null(title)) title else paste("Wing", metric_label)
  p + labs(title=plot_title, subtitle=subtitle, x=lab1, y=lab2) +
    coord_fixed(ratio = 1.05, expand = FALSE) + base_theme
}