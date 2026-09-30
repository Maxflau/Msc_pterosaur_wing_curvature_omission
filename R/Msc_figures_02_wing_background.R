# Build the faint wing-shape background used behind the morphospace plots.
cat("\n--- Building wing profile background ---\n")

# Set how many background wings to draw and how large they appear.
WING_GRID_N <- 10     # profiles per axis; the reference figure uses ~13
WING_SCALE  <- 0.015  # outline size as a fraction of morphospace extent
WING_EXPAND <- 0.05   # same padding as stress_grid

build_wing_lattice <- function() {
  if (!exists("outlines_list") || length(outlines_list) == 0) {
    cat("outlines_list unavailable - no wing background\n"); return(NULL)
  }

  pc1_rng <- range(performance_data_clean$PC1, na.rm = TRUE)
  pc2_rng <- range(performance_data_clean$PC2, na.rm = TRUE)
  p1 <- diff(pc1_rng) * WING_EXPAND
  p2 <- diff(pc2_rng) * WING_EXPAND

  gx <- seq(pc1_rng[1] - p1, pc1_rng[2] + p1, length.out = WING_GRID_N)
  gy <- seq(pc2_rng[1] - p2, pc2_rng[2] + p2, length.out = WING_GRID_N)
  cells <- expand.grid(PC1 = gx, PC2 = gy)

  sc <- WING_SCALE * max(diff(range(gx)), diff(range(gy)))

  out <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    # Nearest specimen supplies the SHAPE; the outline is drawn AT THE CELL.
    # No distance cut-off: every cell gets a profile so the lattice is complete.
    d2 <- (performance_data_clean$PC1 - cells$PC1[i])^2 +
      (performance_data_clean$PC2 - cells$PC2[i])^2
    j <- which.min(d2)
    o <- outlines_list[[performance_data_clean$species[j]]]
    if (is.null(o) || nrow(o) < 3) next

    x <- o$x - mean(o$x); y <- o$y - mean(o$y)
    m <- max(abs(c(x, y)))
    if (!is.finite(m) || m == 0) next

    out[[i]] <- data.frame(PC1 = cells$PC1[i] + sc * x / m,
                           PC2 = cells$PC2[i] + sc * y / m,
                           grid_id = i)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) NULL else do.call(rbind, out)
}

################################################################################
# Choose whether to reuse an existing background or build a new one
#################################################################################
if (exists("wing_overlay_full") && is.data.frame(wing_overlay_full) &&
    nrow(wing_overlay_full) > 0 &&
    all(c("PC1", "PC2", "grid_id") %in% colnames(wing_overlay_full))) {
  wing_bg <- wing_overlay_full
  cat(sprintf("Using wing_overlay_full from Msc_data_visual_v2.R: %d profiles\n",
              length(unique(wing_bg$grid_id))))
} else {
  wing_bg <- build_wing_lattice()
  cat(sprintf("Lattice built here: %d of %d cells (%d x %d)\n",if (is.null(wing_bg)) 0 else length(unique(wing_bg$grid_id)), WING_GRID_N^2, WING_GRID_N, WING_GRID_N))
}

# The reference figure keeps the outlines faint and behind everything else.
# Adjust WING_GRID_N for density and WING_SCALE for size.
if (!is.null(wing_bg)) {
  cat(sprintf("PC1 span %.2f to %.2f | PC2 span %.2f to %.2f\n",
              min(wing_bg$PC1), max(wing_bg$PC1),
              min(wing_bg$PC2), max(wing_bg$PC2)))
}
cat("\n")
