for (dd in c("output/plots", "output/plots_PDF")) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("temporal_data")) stop("temporal_data not found - source Msc_temp_3_1.R.")
if (!"Midpoint" %in% names(temporal_data)) {
  stop("Midpoint absent")
}

# Stacking plot preparations
GROUP_VAR    <- "clade"
PANEL_ASPECT <- 0.45   # plane height / width, in plot units
SHEAR        <- 0.34   # horizontal shift from front to back edge of a plane
PLANE_GAP    <- 0.62   # vertical spacing between planes
HULL_ALPHA   <- 0.38
PT_CEX       <- 0.50

grp_all <- sort(unique(as.character(temporal_data[[GROUP_VAR]])))
grp_all <- grp_all[!is.na(grp_all) & grp_all != ""]
PAL <- if (length(grp_all) <= 12) {
  brewer.pal(max(3, min(12, length(grp_all))), "Set3")[seq_along(grp_all)]
} else colorRampPalette(brewer.pal(12, "Set3"))(length(grp_all))
names(PAL) <- grp_all
cat(sprintf("Groups: %d\n", length(grp_all)))

# ── BUILDER ───────────────────────────────────────────────────────────────────
stacked_panels <- function(d, xv, yv, xlab, ylab, main_title) {
  
  d <- d[is.finite(d[[xv]]) & is.finite(d[[yv]]) & !is.na(d$Time_Bin), ]
  
  #Stratigraphic layout
  bin_age <- tapply(d$Midpoint, d$Time_Bin, mean, na.rm = TRUE)
  bin_age <- bin_age[!is.na(bin_age)]
  bins <- names(sort(bin_age, decreasing = TRUE))
  nb <- length(bins)
  if (nb < 2) { cat("SKIPPED", main_title, "- fewer than two bins\n"); return(invisible(NULL)) }
  
  cat(sprintf("\n  %s\n  Stacking order (bottom to top):\n", main_title))
  for (b in bins) cat(sprintf("    %-24s ~%.0f Ma  (n = %d)\n",
                              b, bin_age[[b]], sum(d$Time_Bin == b)))
  
  # Rescale ONCE on the whole dataset, so every plane shares the same axes
  xr <- range(d[[xv]]); yr <- range(d[[yv]])
  xpad <- diff(xr) * 0.10; ypad <- diff(yr) * 0.10
  xr <- xr + c(-xpad, xpad); yr <- yr + c(-ypad, ypad)
  ux <- (d[[xv]] - xr[1]) / diff(xr)                  # 0 to 1
  uy <- (d[[yv]] - yr[1]) / diff(yr) * PANEL_ASPECT   # 0 to PANEL_ASPECT
  total_h <- PANEL_ASPECT + PLANE_GAP * (nb - 1)
  plot(NA, xlim = c(-0.16, 1 + SHEAR + 0.24), ylim = c(-0.14, total_h + 0.08),
 axes = FALSE, xlab = "", ylab = "", asp = 1)
 title(main_title, line = 0.4, cex.main = 1.2, adj = 0)
  
  for (k in seq_along(bins)) {
    
    base <- PLANE_GAP * (k - 1)
    tx <- function(x, y) x + (y / PANEL_ASPECT) * SHEAR   
    ty <- function(y) base + y
    
    polygon(c(tx(0, 0), tx(1, 0), tx(1, PANEL_ASPECT), tx(0, PANEL_ASPECT)),
 c(ty(0), ty(0), ty(PANEL_ASPECT), ty(PANEL_ASPECT)),
border = "grey25", lwd = 1.1,col = adjustcolor("grey97", alpha.f = 0.85))
    
    sel <- d$Time_Bin == bins[k]
    px <- tx(ux[sel], uy[sel]); py <- ty(uy[sel])
    gv <- as.character(d[[GROUP_VAR]][sel])
    
# Hull preparations
    for (g in unique(gv)) {
      if (is.na(g) || g == "") next
      i <- which(gv == g)
      if (length(i) < 3) next
      h <- chull(px[i], py[i])
      polygon(px[i][h], py[i][h],
              col = adjustcolor(PAL[[g]], alpha.f = HULL_ALPHA),
              border = PAL[[g]], lwd = 1.3)
    }
    
    cols <- PAL[gv]; cols[is.na(cols)] <- "grey70"
    points(px, py, pch = 21, bg = cols, col = "grey15", cex = PT_CEX, lwd = 0.45)
    
    text(-0.035, ty(PANEL_ASPECT / 2),
         labels = sprintf("%s (~%.0f Ma)", bins[k], bin_age[[bins[k]]]),
         srt = 90, cex = 0.68, font = 2, xpd = TRUE)
    text(tx(1, PANEL_ASPECT) - 0.04, ty(PANEL_ASPECT) - 0.035,
         labels = sprintf("n=%d", sum(sel)), cex = 0.62, col = "grey40")
  }
  
  # Ticks on the bottom plane only: all planes share the same axes, so five
  # identical sets would only add clutter
  at_x <- pretty(xr, 5); ax <- (at_x - xr[1]) / diff(xr)
  ok <- ax >= 0 & ax <= 1
  segments(ax[ok], 0, ax[ok], -0.022, xpd = TRUE)
  text(ax[ok], -0.055, labels = at_x[ok], cex = 0.62, xpd = TRUE)
  text(0.5, -0.105, labels = xlab, cex = 0.82, xpd = TRUE)
  
#PC2 axis preparation for fitting part 1
  tx0 <- function(x, y) x + (y / PANEL_ASPECT) * SHEAR
  ty0 <- function(y) y
  
  at_y <- pretty(yr, 4); ay <- (at_y - yr[1]) / diff(yr) * PANEL_ASPECT
  ok <- ay >= 0 & ay <= PANEL_ASPECT
  
#PC2 axis preparation for fitting part 2
  ex <- tx0(1, ay[ok]); ey <- ty0(ay[ok])
  segments(ex, ey, ex + 0.03, ey, xpd = TRUE)
  text(ex + 0.06, ey, labels = at_y[ok], cex = 0.62, xpd = TRUE)
  text(tx0(1, PANEL_ASPECT / 2) + 0.13, PANEL_ASPECT / 2,
       labels = ylab, srt = 90, cex = 0.82, xpd = TRUE)
  
  invisible(TRUE)
}

save_stack <- function(expr, name, w = 9, h = 11.5) {
  for (fmt in c("png", "pdf")) {
    path <- sprintf("output/plots%s/%s.%s",
   if (fmt == "pdf") "_PDF" else "", name, fmt)
    if (fmt == "png") png(path, width = w, height = h, units = "in", res = 300)
    else pdf(path, width = w, height = h)
    layout(matrix(1:2, ncol = 1), heights = c(5.2, 1))
    par(mar = c(1.5, 3.5, 2.5, 1), xpd = TRUE); eval(expr)
    par(mar = c(0, 1, 0, 1)); plot.new()
    legend("center", legend = names(PAL), pt.bg = PAL, pch = 21, col = "grey15",
           ncol = 3, bty = "o", cex = 0.82, pt.cex = 1.3, title = "Clade")
    dev.off()
  }
  cat("  Written:", name, "\n")
}

pct <- function(p, k) if (is.null(p)) "" else
  sprintf(" (%.1f%%)", summary(p)$importance[2, k] * 100)

# ── EFA SHAPE MORPHOSPACE ─────────────────────────────────────────────────────
if (exists("shape_scores") && all(c("shapePC1","shapePC2") %in% names(shape_scores))) {
  sp <- if ("species" %in% names(shape_scores)) "species" else names(shape_scores)[1]
  sd_ <- merge(temporal_data, shape_scores[, c(sp, "shapePC1", "shapePC2")],
               by.x = "species", by.y = sp)
  save_stack(quote(stacked_panels(sd_, "shapePC1", "shapePC2",
                                  paste0("Shape PC1", pct(if (exists("pca_shape")) pca_shape else NULL, 1)),
                                  paste0("Shape PC2", pct(if (exists("pca_shape")) pca_shape else NULL, 2)),
                                  "EFA shape morphospace through time")),
             "STACK_01_shape_morphospace_time")
} else cat("shape_scores unavailable - shape stack skipped.\n")

# ── PERFORMANCE MORPHOSPACE ───────────────────────────────────────────────────
if (all(c("PC1","PC2") %in% names(temporal_data))) {
  save_stack(quote(stacked_panels(temporal_data, "PC1", "PC2",
paste0("PC1", pct(if (exists("pca_performance")) pca_performance else NULL, 1)),
 paste0("PC2", pct(if (exists("pca_performance")) pca_performance else NULL, 2)),
                                  "Performance morphospace through time")),
             "STACK_02_performance_morphospace_time")
} else cat("PC1/PC2 absent - performance stack skipped.\n")