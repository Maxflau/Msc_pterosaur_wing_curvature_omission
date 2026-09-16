# ==============================================================================
# Msc_PGLS.R   —   PHYLOGENETIC TESTS
# ==============================================================================

dir.create("clade_analyses_v4/results", showWarnings = FALSE, recursive = TRUE)

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9_]", "", gsub("[ ]+", "_", trimws(as.character(x)))))
}
if (!exists("performance_data_clean")) stop("performance_data_clean not found.")

# ── 1. TREE ───────────────────────────────────────────────────────────────────
if (!exists("phy_pruned")) {
  cand <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), function(x)
    inherits(get(x, envir = .GlobalEnv), "phylo"), logical(1))]
  if (length(cand) > 0) {
    phy_pruned <- get(cand[1]); cat("Using phylo object:", cand[1], "\n")
  } else if (file.exists("Henry_updated.nex")) {
    tr0 <- read.nexus("Henry_updated.nex")
    phy_pruned <- if (inherits(tr0, "multiPhylo")) tr0[[1]] else tr0
    cat(sprintf("Tree read: %d tips\n", Ntip(phy_pruned)))
  } else stop("No tree available in ", getwd())
}

# ── 2. VARIABLES ──────────────────────────────────────────────────────────────
# Single metric set: von_mises_stress replaces both the old stress_root
# (non-dimensional) and the older pixel-dependent von_mises_stress.
CORE <- c("aspect_ratio", "r2_hat", "von_mises_stress", "wing_loading_ratio")
if (!all(CORE %in% names(performance_data_clean))) {
  stop("Metric set incomplete. Needs: ", paste(CORE, collapse = ", "),
       "\n  Present: ", paste(intersect(CORE, names(performance_data_clean)), collapse = ", "))
}
SET_NAME <- "von_mises_stress"

VARS <- intersect(c(CORE, "wing_curvature","shape_complexity", "pareto_rank_ratio", 
                    "optimality_mean", "optimality_geom", "PC1", "PC2"),
                  colnames(performance_data_clean))
cat("Metric set:", SET_NAME, "\nVariables:", paste(VARS, collapse = ", "), "\n")

# ── 3. ALIGN DATA AND TREE ────────────────────────────────────────────────────
d   <- performance_data_clean
key <- normalise_name(d$species); tip <- normalise_name(phy_pruned$tip.label)
common <- intersect(key, tip)
cat(sprintf("Matched: %d of %d\n\n", length(common), nrow(d)))
if (length(common) < 20) stop("Too few matches.")

d  <- d[key %in% common, ]
tr <- drop.tip(phy_pruned, phy_pruned$tip.label[!tip %in% common])
tr$tip.label <- normalise_name(tr$tip.label)
d  <- d[match(tr$tip.label, normalise_name(d$species)), ]
if (any(tr$edge.length <= 0)) tr$edge.length[tr$edge.length <= 0] <- 1e-6 * max(tr$edge.length)

# ── 4. ROBUST PGLS FIT ────────────────────────────────────────────────────────
#' corPagel optimises from a starting lambda and can fail when the true value is
#' far from it. Starting at 1 failed for pareto_rank_ratio, whose signal lambda
#' is 0.140. Several starts are tried, then Brownian, then OLS.
fit_pgls <- function(df, trs, vname) {
  for (start in c(0.5, 0.9, 0.1, 0.99)) {
    m <- tryCatch(gls(y ~ Midpoint, data = df,
                      correlation = corPagel(start, phy = trs, form = ~1),
                      method = "ML"), error = function(e) NULL)
    if (!is.null(m)) {
      lam <- as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))
      return(list(model = m, lambda = lam, method = sprintf("Pagel (start %.2f)", start)))
    }
  }
  # Brownian fixes lambda at 1 and overcorrects, but it converges where the
  # optimiser will not: reported as such rather than silently substituted
  m <- tryCatch(gls(y ~ Midpoint, data = df,
                    correlation = corBrownian(phy = trs, form = ~1),
                    method = "ML"), error = function(e) NULL)
  if (!is.null(m)) return(list(model = m, lambda = 1, method = "Brownian (fallback)"))
  
  cat("    all PGLS attempts failed for", vname, "\n")
  list(model = NULL, lambda = NA_real_, method = "failed")
}

# ── 5. SIGNAL, OLS AND PGLS ───────────────────────────────────────────────────
res <- do.call(rbind, lapply(VARS, function(v) {
  
  y <- d[[v]]
  
  # Rank ratios and composites are bounded on 0-1, so their residuals cannot be
  # gaussian near the bounds: a logit puts them on an unbounded scale
  if (v %in% c("pareto_rank_ratio", "optimality_mean", "optimality_geom")) {
    y <- qlogis(pmin(pmax(y, 0.001), 0.999))
  } else if (!v %in% c("PC1", "PC2") && all(y > 0, na.rm = TRUE)) {
    y <- log(y)
  }
  
  ok <- is.finite(y) & is.finite(d$Midpoint)
  if (sum(ok) < 20) { cat(sprintf("%-20s skipped (n=%d)\n", v, sum(ok))); return(NULL) }
  
  trs <- drop.tip(tr, tr$tip.label[!ok])
  lam <- tryCatch(phylosig(trs, setNames(y[ok], trs$tip.label),
                           method = "lambda", test = TRUE), error = function(e) NULL)
  
  df <- data.frame(y = y[ok], Midpoint = d$Midpoint[ok])
  o  <- summary(lm(y ~ Midpoint, data = df))$coefficients
  f  <- fit_pgls(df, trs, v)
  
  if (is.null(f$model)) {
    ps <- pse <- pp <- NA_real_
  } else {
    tt <- summary(f$model)$tTable
    ps <- tt["Midpoint","Value"]; pse <- tt["Midpoint","Std.Error"]
    pp <- tt["Midpoint","p-value"]
  }
  
  cat(sprintf("%-20s n=%3d | OLS p=%.4g | PGLS p=%.4g | lambda=%.3f | %s\n",
              v, sum(ok), o[2,4], pp, f$lambda, f$method))
  
  data.frame(variable = v, n = sum(ok),
             transform = if (v %in% c("pareto_rank_ratio","optimality_mean",
                                      "optimality_geom")) "logit"
             else if (v %in% c("PC1","PC2")) "none" else "log",
             signal_lambda = if (!is.null(lam)) round(lam$lambda, 3) else NA,
             signal_p      = if (!is.null(lam)) signif(lam$P, 4) else NA,
             ols_slope = round(o[2,1], 6), ols_p = signif(o[2,4], 4),
             pgls_slope = round(ps, 6), pgls_se = round(pse, 6),
             pgls_p = signif(pp, 4), pgls_lambda = round(f$lambda, 3),
             pgls_method = f$method, stringsAsFactors = FALSE)
}))

if (is.null(res)) stop("No variable produced a usable model.")

# The ratio of slopes measures how much of the uncorrected signal was ancestry
res$slope_ratio <- round(res$ols_slope / res$pgls_slope, 2)
res$verdict <- ifelse(is.na(res$pgls_p), "PGLS failed",
                      ifelse(res$ols_p < 0.05 & res$pgls_p < 0.05, "trend holds",
                             ifelse(res$ols_p < 0.05, "LOST after correction", "no trend either way")))

cat("\n")
print(res[, c("variable","n","transform","signal_lambda","ols_p","pgls_p",
              "pgls_lambda","slope_ratio","verdict")], row.names = FALSE)

res$metric_set <- SET_NAME
write.csv(res, sprintf("clade_analyses_v4/results/pgls_%s_metrics.csv", SET_NAME),row.names = FALSE)
cat(sprintf("\nWritten: pgls_%s_metrics.csv\n", SET_NAME))
cat("Report OLS and PGLS together: the gap between them is the result.\n\n")