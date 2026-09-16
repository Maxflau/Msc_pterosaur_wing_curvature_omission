if (!exists("MV_MATRIX")) stop("Source Msc_stats_00_setup.R first.")
library(vegan)

if (!exists("D_FULL")) D_FULL <- vegdist(MV_MATRIX, method = "euclidean")

N_PERM_PAIR <- 999

# --------------------------------------------------------------------------------
# 1. Pairwise PERMANOVA at both taxonomic scales
# --------------------------------------------------------------------------------
pairwise_permanova <- function(f) {
  
  grp <- as.character(MV_META[[f]])
  ok <- !is.na(grp) & grp != ""
  keep_lv <- names(which(table(grp[ok]) >= MIN_N))
  ok <- ok & grp %in% keep_lv
  lv <- sort(unique(grp[ok]))
  if (length(lv) < 2) return(NULL)
  
  Dm <- as.matrix(D_FULL)
  
  do.call(rbind, lapply(combn(lv, 2, simplify = FALSE), function(pr) {
    sel <- ok & grp %in% pr
    # Six specimens is the practical floor for a two-group permutation test
    if (sum(sel) < 6) return(NULL)
    a <- adonis2(as.dist(Dm[sel, sel]) ~ factor(grp[sel]),
                 permutations = N_PERM_PAIR)
    data.frame(factor = f, group_1 = pr[1], group_2 = pr[2],
               n = sum(sel), R2 = round(a$R2[1], 4),
               F_value = round(a$F[1], 3), p_value = a$`Pr(>F)`[1],
               stringsAsFactors = FALSE)
  }))
}

for (f in intersect(c("clade", "Order"), FACTORS)) {
  
  pw <- pairwise_permanova(f)
  if (is.null(pw)) { cat("SKIPPED", f, "- no usable pair\n"); next }
  
  pw$p_adjusted_BH <- signif(p.adjust(pw$p_value, method = "BH"), 4)
  pw$significant <- pw$p_adjusted_BH < 0.05
  pw <- pw[order(-pw$R2), ]
  
  cat(sprintf("Pairwise PERMANOVA by %s: %d comparisons, %d significant after BH\n",
              f, nrow(pw), sum(pw$significant)))
  write_supp(pw, paste0("S8_pairwise_permanova_by_", tolower(f)))
}

# --------------------------------------------------------------------------------
# 2. PCA loadings — the table that makes the axes interpretable
# --------------------------------------------------------------------------------
# Without this table, PC1 and PC2 are unnamed axes. With it, a reader can see
# which metrics drive each axis and in which direction.
if (exists("pca_performance")) {
  
  imp <- summary(pca_performance)$importance
  n_ax <- min(4, ncol(pca_performance$rotation))
  
  load_tab <- data.frame(
    metric = rownames(pca_performance$rotation),
    round(pca_performance$rotation[, seq_len(n_ax), drop = FALSE], 4),
    stringsAsFactors = FALSE)
  
  # The metric with the largest absolute loading on each axis, for the caption
  dominant <- vapply(seq_len(n_ax), function(k) {
    load_tab$metric[which.max(abs(pca_performance$rotation[, k]))]
  }, character(1))
  
  var_tab <- data.frame(
    axis = colnames(imp)[seq_len(n_ax)],
    proportion_variance = round(imp[2, seq_len(n_ax)], 4),
    cumulative_variance = round(imp[3, seq_len(n_ax)], 4),
    dominant_metric = dominant,
    stringsAsFactors = FALSE)
  
  cat("\nPCA loadings:\n"); print(load_tab, row.names = FALSE)
  cat("\nVariance explained:\n"); print(var_tab, row.names = FALSE)
  
  write_supp(load_tab, "S9_pca_loadings")
  write_supp(var_tab, "S10_pca_variance_explained")
} else {
  cat("\npca_performance not found - loading tables skipped.\n")
}

cat("\nAll supplementary tables written to ", OUT_DIR, "/\n\n", sep = "")