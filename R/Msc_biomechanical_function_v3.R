# --------------------------------------------------------------------------------
# 1. Non-dimensional second moment of area  (r2_hat)
# --------------------------------------------------------------------------------

calculate_r2_raw <- function(outline_coords) {
  x <- outline_coords$x
  y <- outline_coords$y
  n <- length(x)
  j <- c(2:n, 1)
  abs(sum((x * y[j] - x[j] * y) * (x^2 + x * x[j] + x[j]^2)) / 12)
}

anchor_at_base <- function(outline_coords) {
  x <- outline_coords$x
  y <- outline_coords$y
  lo <- min(x); hi <- max(x); R <- hi - lo
  if (!is.finite(R) || R <= 0) return(NULL)
  
  chord_lo <- diff(range(y[x < lo + 0.1 * R]))
  chord_hi <- diff(range(y[x > hi - 0.1 * R]))
  
  if (chord_lo >= chord_hi) {
    data.frame(x = x - lo, y = y)
  } else {
    data.frame(x = hi - x, y = y)
  }
}

calculate_r2_hat <- function(outline_coords) {
  o <- anchor_at_base(outline_coords)
  if (is.null(o)) return(NA_real_)
  
  r2 <- calculate_r2_raw(o)
  S  <- abs(pracma::polyarea(o$x, o$y))
  R  <- diff(range(o$x))
  
  if (!is.finite(S) || S <= 0 || !is.finite(r2) || r2 <= 0) return(NA_real_)
  sqrt(r2 / (S * R^2))
}

# --------------------------------------------------------------------------------
# 2. Von Mises stress (root-section beam model, scale-invariant)
# --------------------------------------------------------------------------------

calculate_von_mises_stress <- function(outline_coords, mass_kg, wingspan_cm) {
  if (is.null(mass_kg) || length(mass_kg) == 0 || is.na(mass_kg) || mass_kg <= 0) {
    mass_kg <- 1.0
  }
  if (is.null(wingspan_cm) || length(wingspan_cm) == 0 || is.na(wingspan_cm) || wingspan_cm <= 0) {
    wingspan_cm <- 50
  }
  
  chord <- max(outline_coords$x) - min(outline_coords$x)
  thickness <- 0.01 * chord
  
  weight_N <- mass_kg * 9.81
  bending_moment <- (weight_N * wingspan_cm / 100) / 4
  
  I <- calculate_second_moment_DIAGNOSTIC(outline_coords)
  
  stress_MPa <- (bending_moment * (thickness/2)) / (I / 10^8)
  
  size_scale <- mass_kg * wingspan_cm
  stress_index <- stress_MPa / size_scale
  
  return(abs(stress_index))
}

# --------------------------------------------------------------------------------
# 3. Aspect ratio computed from the outline itself
# --------------------------------------------------------------------------------

calculate_aspect_ratio <- function(outline_coords) {
  S <- abs(pracma::polyarea(outline_coords$x, outline_coords$y))
  R <- diff(range(outline_coords$x))
  if (!is.finite(S) || S <= 0 || R <= 0) return(NA_real_)
  R^2 / S
}

# --------------------------------------------------------------------------------
# 4. Scale-invariant shape descriptors
# --------------------------------------------------------------------------------

calculate_wing_curvature <- function(outline_coords) {
  le <- which.max(outline_coords$x)
  te <- which.min(outline_coords$x)
  
  x1 <- outline_coords$x[le]; y1 <- outline_coords$y[le]
  x2 <- outline_coords$x[te]; y2 <- outline_coords$y[te]
  
  chord_length <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  if (!is.finite(chord_length) || chord_length == 0) return(NA_real_)
  
  dist <- abs((y2 - y1) * outline_coords$x -
                (x2 - x1) * outline_coords$y +
                x2 * y1 - y2 * x1) / chord_length
  
  max(dist) / chord_length
}

calculate_shape_complexity <- function(outline_coords) {
  dx <- diff(outline_coords$x); dy <- diff(outline_coords$y)
  perimeter <- sum(sqrt(dx^2 + dy^2))
  area <- abs(pracma::polyarea(outline_coords$x, outline_coords$y))
  if (!is.finite(area) || area <= 0) return(NA_real_)
  perimeter^2 / (4 * pi * area)
}

# --------------------------------------------------------------------------------
# 6. Superseded functions — kept for the methods section, never analysed
# --------------------------------------------------------------------------------

calculate_second_moment_DIAGNOSTIC <- function(outline_coords) {
  cx <- mean(outline_coords$x); cy <- mean(outline_coords$y)
  n <- nrow(outline_coords)
  Ix <- 0; Iy <- 0
  for (i in 1:(n - 1)) {
    x1 <- outline_coords$x[i]     - cx; y1 <- outline_coords$y[i]     - cy
    x2 <- outline_coords$x[i + 1] - cx; y2 <- outline_coords$y[i + 1] - cy
    dA <- (x1 * y2 - x2 * y1) / 2
    Ix <- Ix + dA * (y1^2 + y1 * y2 + y2^2) / 3
    Iy <- Iy + dA * (x1^2 + x1 * x2 + x2^2) / 3
  }
  abs(Ix + Iy)
}

cat("Functions defined: calculate_r2_hat(), calculate_von_mises_stress(),\n")