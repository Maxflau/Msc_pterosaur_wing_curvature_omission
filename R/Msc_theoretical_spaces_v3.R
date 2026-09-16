cat("\n================================================================================\n")
cat("PART 5: THEORETICAL MORPHOSPACE FROM EFA COEFFICIENTS\n")
cat("================================================================================\n\n")

if (!exists("pca_shape") || !exists("efa_obj")) {
  stop("pca_shape / efa_obj not found. Source MSc_rep.R before this script.")
}

set.seed(123)
GRID_N      <- 21     # points per axis
EXPAND_FRAC <- 0.20   # extend 20% beyond the empirical range, as in Liu et al.
NB_PTS      <- 100    # points per reconstructed outline


if (exists("outlines_list") && length(outlines_list) > 0) {
  real_spans <- vapply(outlines_list, function(o) {
    if (is.null(o) || nrow(o) < 3) return(NA_real_)
    diff(range(o$x))
  }, numeric(1))
  REF_SPAN <- median(real_spans, na.rm = TRUE)
  cat(sprintf("Reference span for theoretical reconstructions: %.4g (median of %d real outlines)\n",
              REF_SPAN, sum(is.finite(real_spans))))
} else {
  REF_SPAN <- NA_real_
  cat("WARNING: outlines_list not found - theoretical outlines will NOT be\n")
  cat("rescaled to a physical size. von_mises_stress on theoretical shapes\n")
  cat("will not be comparable to real specimens.\n")
}
cat("\n")


####################################################################################
# 0. Pre-flight diagnostics
####################################################################################
rot <- pca_shape$rotation
ctr <- pca_shape$center

if (is.null(rot)) stop("pca_shape has no $rotation — is it a prcomp/PCA object?")
if (is.null(ctr)) { ctr <- rep(0, nrow(rot)); cat("No $center; assuming zero.\n") }

n_coe <- nrow(rot)
nb_h  <- if (!is.null(efa_obj$nb.h)) efa_obj$nb.h else n_coe / 4L

cat(sprintf("Coefficients: %d | harmonics (efa_obj$nb.h): %s\n",
            n_coe, ifelse(is.null(efa_obj$nb.h), "unset", efa_obj$nb.h)))
cat("First coefficient names: ", paste(head(rownames(rot), 8), collapse = ", "), "\n")

# Momocs drops A1/B1/C1/D1 when norm = TRUE, so 4 * nb.h may exceed n_coe
if (n_coe %% 4 != 0) {
  stop("Coefficient count ", n_coe, " is not divisible by 4. Inspect ",
       "colnames(efa_obj$coe) and report it.")
}
nb_h_eff <- n_coe / 4L
if (nb_h_eff != nb_h) {
  cat(sprintf("NOTE: using %d harmonics from the coefficient count, not %s.\n", nb_h_eff, nb_h))
}
nb_h <- nb_h_eff
cat(sprintf("Reconstructing with %d harmonics, %d points per outline\n\n", nb_h, NB_PTS))

####################################################################################
# 1. Grid over the empirical shape-PCA plane
####################################################################################
pc1_rng <- range(pca_shape$x[, 1]); pc2_rng <- range(pca_shape$x[, 2])
pc1_pad <- diff(pc1_rng) * EXPAND_FRAC
pc2_pad <- diff(pc2_rng) * EXPAND_FRAC

grid <- expand.grid(
  shapePC1 = seq(pc1_rng[1] - pc1_pad, pc1_rng[2] + pc1_pad, length.out = GRID_N),
  shapePC2 = seq(pc2_rng[1] - pc2_pad, pc2_rng[2] + pc2_pad, length.out = GRID_N))

cat(sprintf("Grid: %d x %d = %d theoretical shapes\n", GRID_N, GRID_N, nrow(grid)))

####################################################################################
# 2. Inline orientation — no external dependency
####################################################################################
# Centre, rotate onto the principal axis, then put the TIP at max(x) so the base
# sits at min(x), matching the convention the metric functions expect.
orient_inline <- function(o) {
  o$x <- o$x - mean(o$x); o$y <- o$y - mean(o$y)
  pc <- prcomp(cbind(o$x, o$y))
  r  <- as.matrix(cbind(o$x, o$y)) %*% pc$rotation
  o  <- data.frame(x = r[, 1], y = r[, 2])
  
  lo <- min(o$x); hi <- max(o$x); R <- hi - lo
  chord_lo <- diff(range(o$y[o$x < lo + 0.1 * R]))
  chord_hi <- diff(range(o$y[o$x > hi - 0.1 * R]))
  if (chord_hi > chord_lo) o$x <- -o$x        # wider end must be the base
  o
}

# --------------------------------------------------------------------------------
# 3. Reconstruct an outline from a position in shape-PCA space
# --------------------------------------------------------------------------------
split_coefs <- function(coe) {
  if (exists("coeff_split", where = asNamespace("Momocs"), inherits = FALSE)) {
    Momocs::coeff_split(coe, nb.h = nb_h, cph = 4)
  } else {
    b <- split(as.numeric(coe), rep(1:4, each = nb_h))
    list(an = b[[1]], bn = b[[2]], cn = b[[3]], dn = b[[4]])
  }
}

reconstruct_outline <- function(pc1, pc2) {
  scores <- rep(0, ncol(rot)); scores[1] <- pc1; scores[2] <- pc2
  coe <- as.numeric(scores %*% t(rot)) + ctr
  if (anyNA(coe)) stop("NA in reconstructed coefficients")
  
  cf <- split_coefs(coe)
  shp <- Momocs::efourier_i(list(an = cf$an, bn = cf$bn, cn = cf$cn, dn = cf$dn),
                            nb.h = nb_h, nb.pts = NB_PTS)
  o <- data.frame(x = shp[, 1], y = shp[, 2])
  if (nrow(o) < 10 || anyNA(o)) stop("degenerate reconstruction")
  orient_inline(o)
}

# Fail loudly ONCE, before the loop
test_shape <- tryCatch(reconstruct_outline(grid$shapePC1[1], grid$shapePC2[1]), error = function(e)
 stop("Test reconstruction failed: ", e$message, "\nReport this message with the diagnostics above."))
cat(sprintf("Test reconstruction OK: %d points, area = %.4g, span = %.4g\n\n", nrow(test_shape), abs(pracma::polyarea(test_shape$x, test_shape$y)),
            diff(range(test_shape$x))))
