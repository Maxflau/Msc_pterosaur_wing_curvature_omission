# Measure each rebuilt theoretical wing shape and map it onto the Pareto landscape.
# --------------------------------------------------------------------------------
# 4. Flag rebuilt outlines that cross over themselves
# --------------------------------------------------------------------------------
# Liu et al. exclude self-intersecting theoretical outlines as impossible
# designs. This is the only legitimate use of the word "impossible" here.
is_self_intersecting <- function(o) {
  n <- nrow(o); X <- o$x; Y <- o$y
  cr <- function(px, py, qx, qy, rx, ry) (qx-px)*(ry-py) - (qy-py)*(rx-px)
  for (i in 1:(n - 2)) {
    a1x <- X[i]; a1y <- Y[i]; a2x <- X[i+1]; a2y <- Y[i+1]
    for (j in (i + 2):n) {
      if (i == 1 && j == n) next
      k <- if (j == n) 1 else j + 1
      b1x <- X[j]; b1y <- Y[j]; b2x <- X[k]; b2y <- Y[k]
      d1 <- cr(a1x,a1y,a2x,a2y,b1x,b1y); d2 <- cr(a1x,a1y,a2x,a2y,b2x,b2y)
      d3 <- cr(b1x,b1y,b2x,b2y,a1x,a1y); d4 <- cr(b1x,b1y,b2x,b2y,a2x,a2y)
      if (((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0))) return(TRUE)
    }
  }
  FALSE
}

###################################################################################
# 5. Rebuild every theoretical outline and measure its traits
###################################################################################
cat("Reconstructing outlines and computing functional metrics...\n")

theo_list <- vector("list", nrow(grid))
n_impossible <- 0L; n_failed <- 0L; first_error <- NULL
for (k in seq_len(nrow(grid))) {
  o <- tryCatch(reconstruct_outline(grid$shapePC1[k], grid$shapePC2[k]),
                error = function(e) { if (is.null(first_error))
                  first_error <<- conditionMessage(e); NULL })
  if (is.null(o)) { n_failed <- n_failed + 1L; next }
  impossible <- is_self_intersecting(o)
  if (impossible) n_impossible <- n_impossible + 1L
  theo_list[[k]] <- data.frame(
    shapePC1 = grid$shapePC1[k], shapePC2 = grid$shapePC2[k],
    self_intersecting = impossible,
    aspect_ratio     = if (impossible) NA_real_ else calculate_aspect_ratio(o),
    r2_hat           = if (impossible) NA_real_ else calculate_r2_hat(o),
    von_mises_stress = if (impossible) NA_real_ else calculate_von_mises_stress(o, NA_real_, NA_real_),
    wing_curvature   = if (impossible) NA_real_ else calculate_wing_curvature(o),
    shape_complexity = if (impossible) NA_real_ else calculate_shape_complexity(o))

  if (k %% 50 == 0) cat(paste("  ", k, "/", nrow(grid), "\n"))
}

if (!is.null(first_error)) cat("\nFirst reconstruction error:", first_error, "\n")

theoretical_data <- dplyr::bind_rows(theo_list)

if (nrow(theoretical_data) == 0) {
  stop("No theoretical shapes were built. First error: ", ifelse(is.null(first_error), "none recorded", first_error))
}

cat(sprintf("\nTheoretical shapes built: %d\n", nrow(theoretical_data)))
cat(sprintf("Self-intersecting: %d (%.1f%%)\n", n_impossible, 100 * n_impossible / nrow(theoretical_data)))
cat(sprintf("Reconstruction failures: %d\n\n", n_failed))

cat("Theoretical metric ranges (compare with the empirical ones):\n")
print(summary(theoretical_data[, c("aspect_ratio", "r2_hat", "von_mises_stress", "wing_curvature", "shape_complexity")]))
cat("\n")

###################################################################################
# 6. Rank the usable theoretical shapes by Pareto performance
###################################################################################
THEO_OBJECTIVES <- c(aspect_ratio = "max", aspect_ratio = "min",r2_hat = "max", r2_hat = "min",von_mises_stress = "min", von_mises_stress = "max")

goldberg_rank <- function(M) {
  rk <- rep(NA_integer_, nrow(M)); rem <- seq_len(nrow(M)); k <- 0L
  while (length(rem) > 0) {
    sub <- M[rem, , drop = FALSE]
    nd <- vapply(seq_len(nrow(sub)), function(i) {
      ge <- sweep(sub, 2, sub[i, ], ">="); gt <- sweep(sub, 2, sub[i, ], ">")
      !any(apply(ge, 1, all) & apply(gt, 1, any))
    }, logical(1))
    rk[rem[nd]] <- k; rem <- rem[!nd]; k <- k + 1L
  }
  rk
}

usable <- complete.cases(theoretical_data[, names(THEO_OBJECTIVES)]) &
  !theoretical_data$self_intersecting

M <- as.matrix(theoretical_data[usable, names(THEO_OBJECTIVES)])
for (v in names(THEO_OBJECTIVES)) if (THEO_OBJECTIVES[[v]] == "min") M[, v] <- -M[, v]

RO <- goldberg_rank(M); RS <- goldberg_rank(-M)
theoretical_data$pareto_rank_ratio <- NA_real_
theoretical_data$pareto_rank_ratio[usable] <- ifelse(RO + RS == 0, 1, RS / (RO + RS))

cat("Pareto Rank Ratio on theoretical shapes (1 = on the front):\n")
print(summary(theoretical_data$pareto_rank_ratio))
cat(sprintf("Shapes on the front: %d / %d\n\n", sum(RO == 0), length(RO)))

####################################################################################
# 7. Match each real taxon to the nearest point on the theoretical landscape
####################################################################################
if (exists("shape_scores") &&
    all(c("shapePC1", "shapePC2") %in% colnames(shape_scores))) {
  gd <- theoretical_data[usable, ]
  nn <- vapply(seq_len(nrow(shape_scores)), function(i)
    which.min((gd$shapePC1 - shape_scores$shapePC1[i])^2 +
                (gd$shapePC2 - shape_scores$shapePC2[i])^2), integer(1))
  shape_scores$pareto_rank_ratio <- gd$pareto_rank_ratio[nn]

  cat("Pareto Rank Ratio of the real taxa, read off the theoretical landscape:\n")
  print(summary(shape_scores$pareto_rank_ratio))
  cat("\n")
}

dir.create("output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("output/results/supplementals", showWarnings = FALSE, recursive = TRUE)
write.csv(theoretical_data, "output/results/theoretical_shapes.csv",
          row.names = FALSE)

cat("Theoretical morphospace complete.\n\n")
