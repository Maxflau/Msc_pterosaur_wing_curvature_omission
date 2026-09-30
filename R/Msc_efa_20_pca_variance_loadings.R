# EFA shape-space PCA: variance explained + loadings/structure correlations.
# Uses the Momocs PCA object built from efourier() coefficients (pca_shape).
# Mirrors Msc_impossible_regions_v3.R, which does the same two tables for the
# performance-metric PCA, so both morphospaces are reported the same way.

if (!exists("pca_shape")) {
  stop("pca_shape not found - source Msc_rep.R (efourier + PCA) first.")
}

cat(sprintf("Using EFA PCA object 'pca_shape': %d harmonic coefficients, %d specimens.\n",
            nrow(pca_shape$rotation), nrow(pca_shape$x)))

# ── TABLE 1: variance explained by each shape axis ───────────────────────────

efa_pca_summary <- summary(pca_shape)

efa_variance_table <- data.frame(
  Axis           = paste0("shapePC", seq_len(ncol(pca_shape$x))),
  Eigenvalue     = pca_shape$sdev^2,
  Variance_pct   = efa_pca_summary$importance["Proportion of Variance", ] * 100,
  Cumulative_pct = efa_pca_summary$importance["Cumulative Proportion", ] * 100
)

cat("\nVariance explained by each EFA shape axis:\n")

efa_variance_display <- efa_variance_table
efa_variance_display$Eigenvalue     <- format_fixed(efa_variance_table$Eigenvalue, 4)
efa_variance_display$Variance_pct   <- format_fixed(efa_variance_table$Variance_pct, 2)
efa_variance_display$Cumulative_pct <- format_fixed(efa_variance_table$Cumulative_pct, 2)
print(head(efa_variance_display, 10), row.names = FALSE)

n_report <- which(efa_variance_table$Cumulative_pct >= 95)[1]
if (is.na(n_report)) n_report <- ncol(pca_shape$x)

cat(sprintf("\n%d axes needed to reach 95%% of total shape variance (%s%% with first 3).\n",
            n_report,
            format_fixed(efa_variance_table$Cumulative_pct[min(3, nrow(efa_variance_table))], 1)))

save_csv(efa_variance_table, "EFA_pca_variance_table")

# ── TABLE 2: loadings + structure correlations (r, r2) per axis ─────────────
# Rows are Fourier harmonic coefficients (A1, B1, C1, D1, A2, ...), not named
# traits - the harmonics ARE the biologically interpretable unit here, same
# logic as the performance-PCA table (loading * sqrt(eigenvalue) = structure r).

n_axes_report <- 2   # shapePC1 / shapePC2, matching the rest of the pipeline

efa_loadings_raw <- pca_shape$rotation[, seq_len(n_axes_report), drop = FALSE]
efa_eigenvalues  <- pca_shape$sdev[seq_len(n_axes_report)]^2

efa_structure_r  <- sweep(efa_loadings_raw, 2, sqrt(efa_eigenvalues), `*`)
efa_structure_r2 <- efa_structure_r^2

efa_loadings_table <- data.frame(Harmonic_coef = rownames(efa_loadings_raw))

for (i in seq_len(n_axes_report)) {
  efa_loadings_table[[paste0("Loading_shapePC", i)]] <- format_fixed(efa_loadings_raw[, i], 3)
  efa_loadings_table[[paste0("r_shapePC", i)]]        <- format_fixed(efa_structure_r[, i], 3)
  efa_loadings_table[[paste0("r2_shapePC", i)]]       <- format_fixed(efa_structure_r2[, i], 3)
}

# Sort by |r2| on shapePC1 so the harmonics that matter most appear first.
efa_loadings_table <- efa_loadings_table[order(-efa_structure_r2[, 1]), ]

cat("\nTop 10 harmonic coefficients by structure correlation with shapePC1:\n")
print(head(efa_loadings_table, 10), row.names = FALSE)

save_csv(efa_loadings_table, "EFA_pca_loadings_table")

for (i in seq_len(n_axes_report)) {
  top_coef <- rownames(efa_loadings_raw)[which.max(efa_structure_r2[, i])]
  cat(sprintf("shapePC%d: strongest loading is %s (r2 = %s)\n",
              i, top_coef, format_fixed(max(efa_structure_r2[, i]), 3)))
}
cat("\n")
