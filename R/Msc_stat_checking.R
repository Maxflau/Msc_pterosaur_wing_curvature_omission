# Compare each group test with and without the phylogenetic tree correction.
if (!exists("MV_MATRIX")) stop("Source Msc_stats_00_setup.R first.")
if (!exists("phy_pruned")) stop("phy_pruned not found - source Msc_phylo_prep_v3.R first.")
if (!requireNamespace("geomorph", quietly = TRUE)) {
  stop("Package 'geomorph' is required: install.packages(\"geomorph\")")
}
library(geomorph)

N_ITER <- 999

# --------------------------------------------------------------------------------
# 1. Keep only specimens that appear in both the data table and the tree
# --------------------------------------------------------------------------------
key <- normalise_name(MV_META$species)
tip_key <- normalise_name(phy_pruned$tip.label)
common <- intersect(key, tip_key)

cat(sprintf("Specimens with performance data and a tip: %d of %d\n",
            length(common), nrow(MV_MATRIX)))
if (length(common) < 20) stop("Too few matches between the metric matrix and the tree.")

sel <- key %in% common
Y <- MV_MATRIX[sel, , drop = FALSE]
meta <- MV_META[sel, ]
rownames(Y) <- normalise_name(meta$species)

tr <- drop.tip(phy_pruned, phy_pruned$tip.label[!tip_key %in% common])
tr$tip.label <- normalise_name(tr$tip.label)
Y <- Y[tr$tip.label, , drop = FALSE]
meta <- meta[match(tr$tip.label, normalise_name(meta$species)), ]

cat(sprintf("Aligned matrix: %d specimens x %d metrics\n\n", nrow(Y), ncol(Y)))

# --------------------------------------------------------------------------------
# 2. Run the same group test with and without phylogenetic correction
# --------------------------------------------------------------------------------
run_perf <- function(fac, phylogenetic) {
  g <- as.character(meta[[fac]])
  ok <- !is.na(g) & g != ""
  keep_lv <- names(which(table(g[ok]) >= MIN_N))
  ok <- ok & g %in% keep_lv
  if (sum(ok) < 20 || length(unique(g[ok])) < 2) return(NULL)

  Ysub <- Y[ok, , drop = FALSE]
  trs <- drop.tip(tr, tr$tip.label[!ok])
  gdf <- geomorph.data.frame(perf = Ysub, grp = factor(g[ok]), phy = trs)

  fit <- tryCatch({
    if (phylogenetic) {
      procD.pgls(perf ~ grp, phy = phy, data = gdf, iter = N_ITER,
                 print.progress = FALSE)
    } else {
      procD.lm(perf ~ grp, data = gdf, iter = N_ITER, print.progress = FALSE)
    }
  }, error = function(e) { cat("  ", fac, "failed:", e$message, "\n"); NULL })

  if (is.null(fit)) return(NULL)
  a <- fit$aov.table

  data.frame(factor = fac,
             model = if (phylogenetic) "procD.pgls" else "procD.lm (OLS)",
             n = sum(ok), levels = length(unique(g[ok])), df = a$Df[1],
             Rsq = round(a$Rsq[1], 4), F_value = round(a$F[1], 3),
             Z = round(a$Z[1], 3), p_value = a[["Pr(>F)"]][1],
             stringsAsFactors = FALSE)
}

cat("procD.pgls on the performance metrics (Brownian covariance):\n")
pgls_tab <- do.call(rbind, lapply(FACTORS, function(f) {
  r <- run_perf(f, TRUE)
  if (!is.null(r)) cat(sprintf("  %-42s Rsq=%.4f  p=%.4g\n", f, r$Rsq, r$p_value))
  r
}))

cat("\nprocD.lm, no phylogenetic correction, for contrast:\n")
ols_tab <- do.call(rbind, lapply(FACTORS, function(f) {
  r <- run_perf(f, FALSE)
  if (!is.null(r)) cat(sprintf("  %-42s Rsq=%.4f  p=%.4g\n", f, r$Rsq, r$p_value))
  r
}))

if (is.null(pgls_tab)) stop("No factor produced a usable procD.pgls fit.")

pgls_tab$p_adjusted_BH <- signif(p.adjust(pgls_tab$p_value, "BH"), 4)
pgls_tab$significant <- pgls_tab$p_adjusted_BH < 0.05
write_supp(pgls_tab, "S6b_perf_procD_pgls")

# --------------------------------------------------------------------------------
# 3. Compare the corrected and uncorrected results side by side
# --------------------------------------------------------------------------------
# A factor significant under OLS but not under PGLS was a phylogenetic pattern,
# not an ecological one. Stated per row so the write-up cannot quote the
# uncorrected value alone.
if (!is.null(ols_tab)) {
  ols_tab$p_adjusted_BH <- signif(p.adjust(ols_tab$p_value, "BH"), 4)

  cmp <- merge(pgls_tab[, c("factor", "n", "Rsq", "p_adjusted_BH")],
               ols_tab[, c("factor", "Rsq", "p_adjusted_BH")],
               by = "factor", suffixes = c("_pgls", "_ols"))

  cmp$verdict <- ifelse(cmp$p_adjusted_BH_pgls < 0.05,
                        "holds after phylogenetic correction",
                        ifelse(cmp$p_adjusted_BH_ols < 0.05,
                               "lost after correction - phylogenetic, not ecological",
                               "not significant either way"))
  cmp$Rsq_drop <- round(cmp$Rsq_ols - cmp$Rsq_pgls, 4)
  cmp <- cmp[order(-cmp$Rsq_pgls), ]

  cat("\nOLS versus PGLS on the performance matrix:\n")
  print(cmp, row.names = FALSE)
  write_supp(cmp, "S6c_perf_ols_vs_pgls")

  n_lost <- sum(cmp$verdict == "lost after correction - phylogenetic, not ecological")
  if (n_lost > 0) {
    cat(sprintf("\n%d factor(s) lose significance under phylogenetic correction.\n",
                n_lost))
    cat("Those are shared-ancestry patterns and must not be reported as\n")
    cat("ecological effects on aerodynamic performance.\n")
  }
}

cat("\nReport S6b as the primary multivariate test on performance. The PERMANOVA\n")
cat("in S6 assumes independent specimens and is a supplementary approximation.\n\n")
