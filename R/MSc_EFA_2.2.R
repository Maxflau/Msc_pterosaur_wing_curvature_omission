
# ── 2. HELPERS ────────────────────────────────────────────────────────────────
gam_surface <- function(x, y, z, nx = 150, ny = 150) {
  df <- data.frame(x = as.numeric(x), y = as.numeric(y), z = as.numeric(z))
  df <- df[complete.cases(df) & is.finite(df$z), ]
  if (nrow(df) < 10) return(NULL)
  m  <- mgcv::gam(z ~ s(x, y, k = min(60, nrow(df) - 1)), data = df)
  xr <- range(df$x) + diff(range(df$x)) * c(-0.06, 0.06)
  yr <- range(df$y) + diff(range(df$y)) * c(-0.06, 0.06)
  g  <- expand.grid(x = seq(xr[1], xr[2], length.out = nx),
                    y = seq(yr[1], yr[2], length.out = ny))
  g$z <- predict(m, newdata = g)
  g
}

make_hulls <- function(df, group_col) {
  df %>%
    mutate(grp = .data[[group_col]]) %>%
    filter(!is.na(grp)) %>%
    group_by(grp) %>% filter(n() >= 3) %>%
    do({ idx <- chull(.$shapePC1, .$shapePC2); idx <- c(idx, idx[1])
    data.frame(shapePC1 = .$shapePC1[idx], shapePC2 = .$shapePC2[idx]) }) %>%
    ungroup()
}

get_labels <- function(df, n_per_clade = 2) {
  df %>% group_by(clade) %>%
    mutate(dist = sqrt((shapePC1 - mean(shapePC1))^2 +
                         (shapePC2 - mean(shapePC2))^2)) %>%
    slice_max(dist, n = n_per_clade) %>% ungroup()
}

perm_text <- function(pv) {
  sprintf("PERMANOVA: R2 = %s | F = %s | p %s %s", pv$r2, pv$f,
          ifelse(pv$p < 0.001, "<", "="),
          ifelse(pv$p < 0.001, "0.001", sprintf("%.4f", pv$p)))
}

# ── 3. SHARED SURFACES ────────────────────────────────────────────────────────
# Non-finite values are FILTERED, not replaced by 1e-9: a fabricated near-zero
# drags the log surface down wherever a specimen is missing.
ok_stress <- is.finite(shape_perf$von_mises_stress) &
  is.finite(shape_perf$shapePC1) & is.finite(shape_perf$shapePC2)
if (sum(ok_stress) < 20) stop("Only ", sum(ok_stress), " usable specimens.")
cat(sprintf("Stress surface fitted on %d of %d specimens\n",
            sum(ok_stress), nrow(shape_perf)))

stress_bg <- gam_surface(shape_perf$shapePC1[ok_stress],
                         shape_perf$shapePC2[ok_stress],
                         log10(shape_perf$von_mises_stress[ok_stress]))

pareto_bg <- gam_surface(shape_perf$shapePC1, shape_perf$shapePC2,
                         shape_perf$pareto_rank_ratio)

cat("Building EFA grid...\n")
grid_dense <- build_grid_wings(20, 16, 0.25)
if (!exists("impossible_region")) impossible_region <- data.frame()