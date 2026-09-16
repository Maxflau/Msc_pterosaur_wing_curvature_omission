cat("\n================================================================================\n")
cat("PHYLOGENETIC SIGNAL AND PGLS\n")
cat("================================================================================\n\n")

library(nlme)

SIGNAL_VARS <- intersect(
  c("aspect_ratio", "r2_hat", "wing_loading", "von_mises_stress",
    "wing_curvature", "shape_complexity", "PC1", "PC2", "pareto_rank_ratio"),
  colnames(phylo_data))

# ####################################################################################
# 1. Phylogenetic signal
# ####################################################################################
# lambda: 0 = star phylogeny, 1 = Brownian motion.
# K:      <1 = relatives less alike than Brownian expects, >1 = more alike.
signal_results <- data.frame()

for (v in SIGNAL_VARS) {
  x <- phylo_data[[v]]
  # Log the metrics; PC scores and the rank ratio are already on a usable scale
  if (!v %in% c("PC1", "PC2", "pareto_rank_ratio") && all(x > 0, na.rm = TRUE)) {
    x <- log(x)
  }
  names(x) <- phy_pruned$tip.label
  x <- x[is.finite(x)]
  if (length(x) < 10) next
  
  tree_v <- drop.tip(phy_pruned, setdiff(phy_pruned$tip.label, names(x)))
  
  lam <- tryCatch(phylosig(tree_v, x, method = "lambda", test = TRUE),
                  error = function(e) NULL)
  kk  <- tryCatch(phylosig(tree_v, x, method = "K", test = TRUE, nsim = 1000),
                  error = function(e) NULL)
  
  signal_results <- rbind(signal_results, data.frame(
    variable = v, n = length(x),
    lambda   = if (!is.null(lam)) round(lam$lambda, 3) else NA,
    lambda_p = if (!is.null(lam)) signif(lam$P, 3)      else NA,
    K        = if (!is.null(kk))  round(kk$K, 3)        else NA,
    K_p      = if (!is.null(kk))  signif(kk$P, 3)       else NA,
    stringsAsFactors = FALSE))
}

cat("PHYLOGENETIC SIGNAL:\n")
print(signal_results, row.names = FALSE)
cat("\nlambda near 1 with a significant p means related taxa resemble each other.\n")
cat("K below 1 means evolution is more labile than Brownian motion, which is\n")
cat("compatible with ecological convergence on a phylogenetically structured\n")
cat("background. Report both, they answer different questions.\n\n")

# --------------------------------------------------------------------------------
# 2. PGLS — does any temporal trend survive phylogenetic correction?
# --------------------------------------------------------------------------------
run_pgls <- function(response) {
  if (!all(c(response, "Midpoint") %in% colnames(phylo_data))) return(invisible(NULL))
  
  df <- phylo_data[, c(response, "Midpoint")]
  names(df)[1] <- "y"
  ok <- complete.cases(df)
  df <- df[ok, ]
  tr <- drop.tip(phy_pruned, phy_pruned$tip.label[!ok])
  
  cat(sprintf("--- %s ~ Midpoint ---\n", response))
  cat("Uncorrected OLS:\n")
  print(round(summary(lm(y ~ Midpoint, data = df))$coefficients, 5))
  
  m_lam <- tryCatch(
    gls(y ~ Midpoint, data = df,
        correlation = corPagel(1, phy = tr, form = ~1), method = "ML"),
    error = function(e) { cat("corPagel failed:", e$message, "\n"); NULL })
  
  if (!is.null(m_lam)) {
    cat("\nPGLS, Pagel's lambda estimated:\n")
    print(round(summary(m_lam)$tTable, 5))
    cat(sprintf("Estimated lambda: %.3f\n\n",
                as.numeric(coef(m_lam$modelStruct$corStruct, unconstrained = FALSE))))
  }
  invisible(m_lam)
}

for (resp in intersect(c("pareto_rank_ratio", "PC1", "aspect_ratio","PC2", "r2_hat","von_mises_stress", "wing_loading_ratio"),
                       colnames(phylo_data))) {
  run_pgls(resp)
}


# ####################################################################################
# 3. Clade differences against the phylogenetic null
# ####################################################################################

# An ordinary ANOVA or Fisher test on clades is partly circular: clades are
# defined on the tree. phylANOVA compares against the correct null.
if ("clade" %in% colnames(phylo_data)) {
  y   <- setNames(phylo_data$PC1, phy_pruned$tip.label)
  grp <- setNames(as.factor(phylo_data$clade), phy_pruned$tip.label)
  ok  <- is.finite(y) & !is.na(grp)
  
  pa <- tryCatch(
    phylANOVA(drop.tip(phy_pruned, phy_pruned$tip.label[!ok]),
              grp[ok], y[ok], nsim = 1000, posthoc = FALSE),
    error = function(e) { cat("phylANOVA failed:", e$message, "\n"); NULL })
  
  if (!is.null(pa)) {
    cat(sprintf("PHYLOGENETIC ANOVA (PC1 ~ clade): F = %.3f, p = %.4f\n\n",
                pa$F, pa$Pf))
  }
}

# ####################################################################################
# 4. Robustness across grafted trees
# ####################################################################################
if (!is.null(all_trees) && length(all_trees) > 1 &&
    all(c("pareto_rank_ratio", "Midpoint") %in% colnames(phylo_data))) {
  
  n_trees <- min(length(all_trees), 100)
  cat(sprintf("Repeating the PGLS across %d trees...\n", n_trees))
  
  slopes <- numeric(0); pvals <- numeric(0)
  for (i in seq_len(n_trees)) {
    tr <- all_trees[[i]]
    kt <- normalise_name(tr$tip.label) %in% phylo_data$tip_key
    tr <- drop.tip(tr, tr$tip.label[!kt])
    dd <- phylo_data[normalise_name(tr$tip.label), c("pareto_rank_ratio", "Midpoint")]
    names(dd)[1] <- "y"
    if (any(!complete.cases(dd))) next
    m <- tryCatch(gls(y ~ Midpoint, data = dd,
                      correlation = corPagel(1, phy = tr, form = ~1), method = "ML"),
                  error = function(e) NULL)
    if (!is.null(m)) {
      slopes <- c(slopes, summary(m)$tTable["Midpoint", "Value"])
      pvals  <- c(pvals,  summary(m)$tTable["Midpoint", "p-value"])
    }
  }
  
  if (length(slopes) > 0) {
    cat(sprintf("Successful fits: %d\n", length(slopes)))
    print(summary(slopes))
    cat(sprintf("Trees with p < 0.05: %d / %d (%.1f%%)\n\n",
                sum(pvals < 0.05), length(pvals), 100 * mean(pvals < 0.05)))
    cat("Report this proportion: it is the honest measure of how much the result\n")
    cat("depends on the placement of the 25 grafted taxa.\n\n")
  }
}

phylo_results <- list(signal = signal_results, data = phylo_data, tree = phy_pruned)
cat("Phylogenetic analysis complete.\n\n")