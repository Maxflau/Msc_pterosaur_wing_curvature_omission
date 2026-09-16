# Build the performance morphospace and mark the observed envelope of points.
cat("\n================================================================================\n")
cat("PART 4: PERFORMANCE MORPHOSPACE AND EMPIRICAL ENVELOPE\n")
cat("================================================================================\n\n")

################################################################################################
# Pick the measurements used to build the morphospace
################################################################################################
CORE_VARS <- c("aspect_ratio", "r2_hat","von_mises_stress", "wing_curvature", "shape_complexity",
               "wing_loading_ratio")

missing_vars <- setdiff(CORE_VARS, colnames(performance_data))
if (length(missing_vars) > 0) {
  stop("Missing performance columns: ", paste(missing_vars, collapse = ", "),
       "\n  Present: ", paste(intersect(
         c("aspect_ratio","r2_hat","von_mises_stress","wing_loading_ratio",
           "wing_curvature","shape_complexity"),
         colnames(performance_data)), collapse = ", "),
       "\n  Run Msc_performance_metrics_v3.R first.")
}

PERFORMANCE_VARS <- CORE_VARS

# ── 2. Build the measurement matrix used for the PCA ──────────────────────────
performance_matrix <- performance_data %>%
  dplyr::select(all_of(PERFORMANCE_VARS)) %>%
  as.matrix()

complete_rows <- complete.cases(performance_matrix)
performance_matrix_clean <- performance_matrix[complete_rows, ]
performance_data_clean   <- performance_data[complete_rows, ]

cat(paste0("Specimens with complete performance data: ",
           nrow(performance_matrix_clean), " / ", nrow(performance_data), "\n\n"))

if (any(performance_matrix_clean <= 0, na.rm = TRUE)) {
  bad <- colnames(performance_matrix_clean)[
    apply(performance_matrix_clean <= 0, 2, any)]
  stop("Non-positive values in: ", paste(bad, collapse = ", "),
       "\n  A zero usually means a degenerate outline.")
}
performance_matrix_log <- log(performance_matrix_clean)

# ── 3. Check whether any measurements are almost duplicates ───────────────────
cm <- cor(performance_matrix_log, method = "spearman")
high <- which(abs(cm) > 0.9 & upper.tri(cm), arr.ind = TRUE)
if (nrow(high) > 0) {
  cat("WARNING - pairs correlating above 0.9:\n")
  for (i in seq_len(nrow(high))) {
    cat(sprintf("  %s / %s : rho = %+.3f\n", rownames(cm)[high[i,1]],
                colnames(cm)[high[i,2]], cm[high[i,1], high[i,2]]))
  }
  cat("Consider dropping one of each pair before interpreting the axes.\n\n")
}

# ── 4. Run the PCA and store the first two axes ───────────────────────────────
pca_performance <- prcomp(performance_matrix_log, scale. = TRUE, center = TRUE)

cat("Performance morphospace PCA:\n")
print(summary(pca_performance)$importance[, 1:min(4, ncol(pca_performance$x))])
cat("\n")

performance_data_clean$PC1 <- pca_performance$x[, 1]
performance_data_clean$PC2 <- pca_performance$x[, 2]

cat("PC1 / PC2 loadings - report these so the axes are interpretable:\n")
print(round(pca_performance$rotation[, 1:2], 3))
cat("\n")

cat("Correlation among log-transformed performance variables:\n")
print(round(cor(performance_matrix_log), 2))
cat("\n")

# ── 5. Build the envelope that contains most observed specimens ───────────────
ENVELOPE_LEVEL <- 0.95

pc_plane <- cbind(performance_data_clean$PC1, performance_data_clean$PC2)
envelope_center <- colMeans(pc_plane)
envelope_cov    <- cov(pc_plane)
envelope_cut    <- qchisq(ENVELOPE_LEVEL, df = 2)

performance_data_clean$mahalanobis_d2 <-
  mahalanobis(pc_plane, center = envelope_center, cov = envelope_cov)

performance_data_clean$morphometric_outlier <-
  performance_data_clean$mahalanobis_d2 > envelope_cut

cat(sprintf("Empirical envelope: %.0f%% chi-square ellipse (cutoff %.2f)\n",
            ENVELOPE_LEVEL * 100, envelope_cut))
cat(sprintf("Morphometric outliers: %d / %d (%.1f%%)\n\n",
            sum(performance_data_clean$morphometric_outlier),
            nrow(performance_data_clean),
            100 * mean(performance_data_clean$morphometric_outlier)))

if (any(performance_data_clean$morphometric_outlier)) {
  cat("Specimens outside the envelope (inspect silhouettes and CSV rows):\n")
  print(performance_data_clean %>%
          dplyr::filter(morphometric_outlier) %>%
          dplyr::arrange(desc(mahalanobis_d2)) %>%
          dplyr::select(species, dplyr::any_of("clade"),
                        mahalanobis_d2, PC1, PC2) %>%
          head(15))
  cat("\n")
}

# ── 6. Save the ellipse shape and the inside/outside test ---------------------
envelope_ellipse <- local({
  theta <- seq(0, 2 * pi, length.out = 200)
  circle <- cbind(cos(theta), sin(theta)) * sqrt(envelope_cut)
  ev <- eigen(envelope_cov)
  pts <- circle %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
  data.frame(PC1 = pts[, 1] + envelope_center[1],
             PC2 = pts[, 2] + envelope_center[2])
})

inside_envelope <- function(pc1, pc2) {
  d2 <- mahalanobis(cbind(pc1, pc2), center = envelope_center, cov = envelope_cov)
  d2 <= envelope_cut
}

cat("Empirical envelope defined. inside_envelope() and envelope_ellipse available.\n\n")
