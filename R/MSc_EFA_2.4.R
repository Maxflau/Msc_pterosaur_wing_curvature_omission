metrics <- list(
  list(col="aspect_ratio",      label="Aspect Ratio\n(WAR)",
       log=FALSE, pal="plasma",  rev=FALSE, brk=c(6,7.5,10,12.5,15),
       clade_col=TRUE,  transp=FALSE, dark_bg=FALSE),
  list(col="r2_hat",            label="Second moment\nof area",
       log=FALSE, pal="viridis", rev=FALSE, brk=NULL,
       clade_col=TRUE,  transp=FALSE, dark_bg=FALSE),
  list(col="von_mises_stress",  label="Von Mises stress\n(log10)",
       log=TRUE,  pal="magma",   rev=FALSE, brk=NULL,
       clade_col=TRUE,  transp=TRUE,  dark_bg=TRUE),
  list(col="wing_curvature",    label="Wing curvature",
       log=FALSE, pal="cividis", rev=TRUE,  brk=NULL,
       clade_col=FALSE, transp=FALSE, dark_bg=FALSE),
  list(col="pareto_rank_ratio", label="Pareto rank ratio",
       log=FALSE, pal="Greens",  rev=FALSE, brk=NULL,
       clade_col=FALSE, transp=FALSE, dark_bg=FALSE))

# ── 7. PRODUCE FIGURES — one plot per data type, pipeline naming ──────────────
dir.create("clade_analyses_v4/plots_PDF", showWarnings=FALSE, recursive=TRUE)
dir.create("clade_analyses_v4/plots",     showWarnings=FALSE, recursive=TRUE)

# Helper: save both PDF and PNG
save_plot <- function(p, name, w=11, h=8.5) {
  base <- file.path("clade_analyses_v4/plots_PDF", name)
  ggsave(paste0(base,".pdf"), p, width=w, height=h, dpi=300)
  ggsave(file.path("clade_analyses_v4/plots", paste0(name,".png")),
         p, width=w, height=h, dpi=150)
  cat(sprintf("✓ Plot saved: %s (PNG + PDF)\n", name))
}

# ── 7A. METRIC MORPHOSPACE PLOTS — 1A / 1B / 1C  ─────────────────────────────
# WAR (1A) and 2MA (1B): clade border + metric fill
# SVM (1C): dark navy bg + transparent points + clade borders
# Others (1D-1E): black hulls, grey point border
metric_labels <- c("1A","1B","1C","1D","1E")

for (i in seq_along(metrics)) {
  m    <- metrics[[i]]
  lbl  <- metric_labels[i]
  p <- make_metric_plot(
    metric_col      = m$col,
    metric_label    = m$label,
    metric_log      = m$log,
    metric_pal      = m$pal,
    metric_rev      = m$rev,
    metric_breaks   = m$brk,
    group_var       = "clade",
    group_cols      = CLADE_COLS,
    hull_fill       = FALSE,
    clade_colour    = m$clade_col,
    transparent_pts = m$transp,
    dark_bg         = m$dark_bg,
    legend_bottom   = FALSE,
    perm_factor     = "clade",
    title           = sprintf("Plot %s: %s with smaller filled wing shapes", lbl, m$label)
  )
  if (is.null(p)) next
  save_plot(p, sprintf("CLADE_0%s_%s", lbl, m$col))
}

# ── 7B. GROUP MORPHOSPACE PLOTS — CLADE_10 to CLADE_14 ───────────────────────

make_group_morph <- function(group_var, group_cols, title, subtitle="",
                             hull_alpha=0.18, legend_bottom=FALSE) {
  hulls <- make_hulls(shape_perf, group_var)
  pts   <- shape_perf %>%
    mutate(grp = .data[[group_var]]) %>%
    filter(!is.na(grp))
  legend_pos <- if (legend_bottom) "bottom" else "right"
  
  p <- ggplot()
  
  # 1. Stress background (inverted Blues: light = high stress)
  if (!is.null(stress_bg)) {
    p <- p +
      geom_raster(data=stress_bg, aes(x,y,fill=z), interpolate=TRUE) +
      scale_fill_distiller(palette = "Blues", direction = 1,
                           name = "Von Mises stress\n(log10)",
                           labels = scales::label_number(accuracy = 0.01),
                           guide = guide_colorbar(order = 99, barheight = unit(3.5,"cm"),
                                                  barwidth = unit(0.32,"cm"),
                                                  title.theme = element_text(size = 7, face = "bold")))
  }
  # 2. Impossible region
  if (nrow(impossible_region) > 0)
    p <- p + geom_tile(data=impossible_region, aes(x,y), fill="grey80", alpha=0.75)
  
  # 3. EFA grid — smaller (0.25) and filled white
  p <- p +
    new_scale_fill() +
    geom_polygon(data=grid_dense, aes(gx,gy,group=grid_id),
                 fill="white", colour="grey45", linewidth=0.06, alpha=0.50)
  
  # 4. Convex hulls: semi-transparent coloured fill + matching border
  p <- p +
    geom_polygon(data=hulls,
                 aes(shapePC1, shapePC2, group=grp, fill=grp, colour=grp),
                 alpha=hull_alpha, linewidth=0.60, show.legend=FALSE) +
    scale_fill_manual(values=group_cols, guide="none") +
    scale_colour_manual(values=group_cols, guide="none")
  
  # 5. Specimen points: WHITE interior, group-coloured BORDER
  #    Matches reference: open circles with category-coloured outline
  #    Legend: same circles with coloured borders, white fill
  p <- p +
    new_scale_colour() +
    geom_point(data=pts,
               aes(shapePC1, shapePC2, colour=grp),
               shape=21, fill=NA,         # transparent interior
               size=2.5, stroke=0.80) +
    scale_colour_manual(values=group_cols,
                        name=gsub("_", " ", group_var),
                        guide=guide_legend(order=1,
                                           ncol=if(legend_bottom) 3L else 1L,
                                           override.aes=list(shape=21,
                                                             fill=NA,
                                                             size=4.0,
                                                             stroke=1.0)))
  
  p + labs(title=title, subtitle=subtitle, x=lab1, y=lab2) +
    coord_fixed(ratio = 0.95, expand = FALSE) + base_theme  # legend always right
}

cat("\n── Group morphospace plots (CLADE_10–14)\n")

save_plot(
  make_group_morph("clade", CLADE_COLS,
                   "Morphospace Occupation by Clade",
                   "Convex hulls by clade"),
  "CLADE_10_morphospace_clade_convex_hulls"
)