# Ensure res_dir exists (may not if Section 9 was skipped)
if (!exists("res_dir")) {
  res_dir <- "output/results/Supplemental_EFA"
  dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
}

required_time_cols <- c("Time_Bin", "Midpoint")
missing_time <- required_time_cols[!required_time_cols %in% names(shape_perf)]

if (length(missing_time) > 0) {
  cat(sprintf("  Missing: %s — skipping Section 10.\n",
              paste(missing_time, collapse=", ")))
} else {
  # ── 10B. PREPARE DATA ─────────────────────────────────────────────────────────
# Keep only the rows needed for the through-time figures.
  bin_order <- c("Late Triassic", "Early+Middle Jurassic", "Late Jurassic",
                 "Early Cretaceous", "Late Cretaceous")
  bin_order <- bin_order[bin_order %in% unique(shape_perf$Time_Bin)]

  shape_temporal <- shape_perf %>%
    filter(!is.na(Time_Bin), !is.na(shapePC1), !is.na(shapePC2), !is.na(clade)) %>%
    mutate(Time_Bin = factor(Time_Bin, levels=bin_order))

  n_bins_efa <- nlevels(droplevels(shape_temporal$Time_Bin))
  cat(sprintf("  %d specimens × %d bins: %s\n",
              nrow(shape_temporal), n_bins_efa,
              paste(levels(droplevels(shape_temporal$Time_Bin)), collapse=" | ")))

  # ── 10C. HULLS ────────────────────────────────────────────────────────────────
# Draw a hull around each clade inside each time bin.
  efa_hulls_time <- shape_temporal %>%
    group_by(Time_Bin, clade) %>%
    filter(n() >= 3) %>%
    do({
      idx <- chull(.$shapePC1, .$shapePC2); idx <- c(idx, idx[1])
      data.frame(shapePC1=.$shapePC1[idx], shapePC2=.$shapePC2[idx])
    }) %>% ungroup()

  # ── 10D. COLOURS ──────────────────────────────────────────────────────────────
# Pick colours for all clades shown through time.
  clades_in_time <- sort(unique(shape_temporal$clade))
  efa_time_cols  <- CLADE_COLS[clades_in_time]
  missing_idx    <- is.na(efa_time_cols)
  if (any(missing_idx)) {
    extra <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(sum(missing_idx))
    efa_time_cols[missing_idx] <- extra
  }

  # ── 10E. PANEL LABELS (time bin + n specimens + n clades) ─────────────────────
# Build clear panel labels with age and sample counts.
  bin_labels <- shape_temporal %>%
    group_by(Time_Bin) %>%
    summarise(mid  = round(mean(Midpoint, na.rm=TRUE), 0),
              n    = n(),
              ncl  = n_distinct(clade), .groups="drop") %>%
    mutate(label = paste0(as.character(Time_Bin),
                          "\n(~", mid, " Ma, n=", n,
                          ", clades=", ncl, ")"))

  shape_temporal <- shape_temporal %>%
    left_join(bin_labels %>% select(Time_Bin, label), by="Time_Bin") %>%
    mutate(Time_label = factor(label, levels=bin_labels$label))

  efa_hulls_time <- efa_hulls_time %>%
    left_join(bin_labels %>% select(Time_Bin, label), by="Time_Bin") %>%
    mutate(Time_label = factor(label, levels=bin_labels$label))

  # ── 10F. GGPLOT — matching reference style ────────────────────────────────────
# Draw the faceted plot that shows shape space through time.
  # White background, grey backdrop, semi-transparent hulls, filled solid circles
  efa_time_plot <- ggplot() +
    # 1. All specimens as small grey backdrop dots
    geom_point(data = shape_perf %>% filter(!is.na(shapePC1)), aes(shapePC1, shapePC2), colour="grey82", size=0.6, alpha=0.55, shape=16) +
    # 2. Semi-transparent filled convex hulls (matching reference alpha=0.3)
    geom_polygon(data = efa_hulls_time, aes(shapePC1, shapePC2, group=interaction(Time_label, clade), fill=clade, colour=clade), alpha=0.28, linewidth=0.65) +
    # 3. Filled solid circles: fill = clade colour, border = black (base R pch=21 style)
    geom_point(data = shape_temporal,
               aes(shapePC1, shapePC2, fill=clade),
               shape=21, colour="black", size=2.2, stroke=0.45) +
    # Scales
    scale_fill_manual(values=efa_time_cols,   name="Clade",
                      guide=guide_legend(ncol=2,
                                         override.aes=list(shape=21, colour="black",
                                                           size=3.5, stroke=0.6))) +
    scale_colour_manual(values=efa_time_cols, guide="none") +
    # Facets: 3 columns, matching reference grid
    facet_wrap(~ Time_label, ncol=3) +
    labs(title    = "EFA Morphospace Through Geological Time",
         subtitle = sprintf("%d specimens across %d time bins",
                            nrow(shape_temporal), n_bins_efa),
         x = lab1, y = lab2) +
    coord_equal() +
    theme_classic(base_size=10) +
    theme(
      strip.text       = element_text(face="bold", size=8.5),
      strip.background = element_rect(fill="grey95", colour="grey60", linewidth=0.4),
      legend.position  = "right",
      legend.text      = element_text(size=7.5),
      legend.title     = element_text(size=9, face="bold"),
      legend.key.size  = unit(0.4, "cm"),
      panel.border     = element_rect(colour="grey40", fill=NA, linewidth=0.4),
      panel.grid.major = element_line(colour="grey90", linewidth=0.2),
      panel.spacing    = unit(0.5, "cm"),
      plot.background  = element_rect(fill="white", colour=NA)
    )

  # ── 10G. DISPARITY THROUGH TIME ───────────────────────────────────────────────
# Summarise how spread-out the shapes are in each time bin.
  efa_disparity <- shape_temporal %>%
    group_by(Time_Bin) %>%
    summarise(n        = n(),
              n_clades = n_distinct(clade),
              SoV_PC1  = var(shapePC1,        na.rm=TRUE),
              SoV_PC2  = var(shapePC2,        na.rm=TRUE),
              disparity= SoV_PC1 + SoV_PC2,
              mean_opt = mean(pareto_rank_ratio, na.rm=TRUE),
              Midpoint = mean(Midpoint,        na.rm=TRUE),
              .groups  = "drop") %>%
    arrange(Midpoint)

  disp_time_plot <- ggplot(efa_disparity, aes(x=Midpoint, y=disparity)) +
    geom_line(colour="#2980B9", linewidth=0.9) +
    geom_point(aes(size=n, fill=mean_opt), shape=21, colour="grey20", stroke=0.5) +
    geom_text(aes(label=paste0(as.character(Time_Bin),"\n(n=",n,")")),
              vjust=-0.7, size=2.4, colour="grey30") +
    scale_fill_gradient(low="white", high="#1A5C2A",
                        name="Mean Pareto\nOptimality", limits=c(0,1)) +
    scale_size_continuous(name="N specimens", range=c(3,10)) +
    scale_x_reverse(name="Age (Ma)") +
    labs(title    = "EFA Morphospace Disparity Through Time",
         subtitle = "Sum of variances (shapePC1 + shapePC2)",
         y = "Disparity (SoV)") +
    theme_bw(base_size=11)

  # ── 10H. SAVE ─────────────────────────────────────────────────────────────────
# Save the plots and the summary table.
  n_rows_efa <- ceiling(n_bins_efa / 3)
  fig_h_efa  <- max(8, n_rows_efa * 5.2)

  dir.create("output/plots",     showWarnings=FALSE, recursive=TRUE)
  dir.create("output/plots_PDF", showWarnings=FALSE, recursive=TRUE)

  ggsave("output/plots_PDF/CLADE_EFA_morphospace_time.pdf",
         efa_time_plot, width=18, height=fig_h_efa, dpi=300)

  ggsave("output/plots_PDF/CLADE_EFA_disparity_time.pdf",
         disp_time_plot, width=10, height=6, dpi=300)

  write.csv(efa_disparity,file.path(res_dir, "efa_disparity_through_time.csv"), row.names=FALSE)
  cat("  ✓ efa_disparity_through_time.csv\n")
} # end if required columns present
