cat("\n=== PARETO FRONT ===\n\n")
# Make sure the output folder exists and the input table is loaded.
dir.create("supplementals", showWarnings = FALSE, recursive = TRUE)
if (!exists("performance_data_clean")) stop("performance_data_clean not found.")

d <- performance_data_clean

# ── 1. OBJECTIVES ──────────────────────────────────────────────────────────────
# List the measures used to judge performance.
# Each one must appear once, with the direction kept the same throughout.
OBJECTIVES <- c(aspect_ratio = "max",aspect_ratio = "min", r2_hat = "max",r2_hat = "min",
                von_mises_stress = "max",von_mises_stress = "min",wing_loading_ratio ="min")

have_obj <- all(names(OBJECTIVES) %in% names(d))

if (!have_obj && all(c("aspect_ratio","r2_hat") %in% names(d)) &&
    exists("calculate_von_mises_stress") && exists("outlines_list")) {
  d$von_mises_stress <- vapply(seq_len(nrow(d)), function(i) {
    s <- d$species[i]
    o <- outlines_list[[s]]
    if (is.null(o) || nrow(o) < 4) return(NA_real_)
    mk <- if ("Mass_kg" %in% colnames(d)) d$Mass_kg[i] else NA_real_
    ws <- if ("Wingspan_cm" %in% colnames(d)) d$Wingspan_cm[i] else NA_real_
    tryCatch(calculate_von_mises_stress(o, mk, ws), error = function(e) NA_real_)
  }, numeric(1))
  have_obj <- sum(is.finite(d$von_mises_stress)) > 20
  if (have_obj) cat(sprintf("von_mises_stress computed here: %d specimens\n",
                            sum(is.finite(d$von_mises_stress))))
}

if (!have_obj) {
  stop("Metric set incomplete.\n",
       "  Present: ", paste(intersect(names(OBJECTIVES), names(d)), collapse = ", "),
       "\n  Needs   : ", paste(names(OBJECTIVES), collapse = ", "))
}
cat("Metric set: aspect_ratio, r2_hat, von_mises_stress, wing_loading_ratio\n")

# ── 2. BUILD THE OBJECTIVE MATRIX ─────────────────────────────────────────────
# Keep only rows that have usable values for every chosen measure.
M <- as.matrix(d[, names(OBJECTIVES)])
keep <- complete.cases(M) & apply(is.finite(M), 1, all)
dd <- d[keep, ]; M <- M[keep, , drop = FALSE]

for (v in names(OBJECTIVES)) if (OBJECTIVES[[v]] == "min") M[, v] <- -M[, v]

cat(sprintf("Specimens: %d of %d | objectives: %s\n\n", nrow(M), nrow(d),
            paste(names(OBJECTIVES), OBJECTIVES, sep = "=", collapse = ", ")))

# ── 3. GOLDBERG RANKS AND PARETO RANK RATIO (Deakin et al. 2022) ──────────────
# Find the best trade-off set, then score how close each specimen is to it.
is_nd <- function(m) vapply(seq_len(nrow(m)), function(i) {
  ge <- sweep(m, 2, m[i, ], ">="); gt <- sweep(m, 2, m[i, ], ">")
  !any(apply(ge, 1, all) & apply(gt, 1, any))
}, logical(1))

goldberg <- function(m) {
  rk <- rep(NA_integer_, nrow(m)); rem <- seq_len(nrow(m)); k <- 0L
  while (length(rem) > 0) {
    nd <- is_nd(m[rem, , drop = FALSE]); rk[rem[nd]] <- k
    rem <- rem[!nd]; k <- k + 1L
  }
  rk
}

RO <- goldberg(M); RS <- goldberg(-M)
dd$goldberg_rank     <- RO
dd$pareto_rank_ratio <- ifelse(RO + RS == 0, 1, RS / (RO + RS))
dd$on_front          <- RO == 0
dd$dominance_count   <- vapply(seq_len(nrow(M)), function(i) {
  ge <- sweep(M, 2, M[i, ], ">="); gt <- sweep(M, 2, M[i, ], ">")
  sum(apply(ge, 1, all) & apply(gt, 1, any))
}, numeric(1))

cat(sprintf("Pareto front: %d / %d (%.1f%%) | ranks 0 to %d\n",
            sum(dd$on_front), nrow(dd), 100 * mean(dd$on_front), max(RO)))

# ── 4. LIMITING OBJECTIVE AND STRATEGY ────────────────────────────────────────
# Work out which measure holds each specimen back the most.
pct <- apply(M, 2, function(x) rank(x, ties.method = "average") / length(x))
colnames(pct) <- names(OBJECTIVES)

dd$limiting_objective  <- colnames(pct)[apply(pct, 1, which.min)]
dd$limiting_percentile <- apply(pct, 1, min)
dd$balance             <- apply(pct, 1, function(p) 1 - (max(p) - min(p)))

for (v in colnames(pct)) dd[[paste0("pareto_score_", v)]] <- pct[, v]

TRADEOFF_CORE <- names(OBJECTIVES)
for (v in TRADEOFF_CORE) {
  others <- setdiff(TRADEOFF_CORE, v)
  dd[[paste0("pareto_tradeoff_", v)]] <- pct[, v] - rowMeans(pct[, others, drop = FALSE])
}

cat("\nPercentile correlations among the trade-off objectives:\n")
print(round(cor(pct[, TRADEOFF_CORE], method = "spearman"), 3))

Ms <- scale(M); front_s <- Ms[dd$on_front, , drop = FALSE]
dd$distance_to_front <- vapply(seq_len(nrow(Ms)), function(i)
  min(sqrt(rowSums((sweep(front_s, 2, Ms[i, ], "-"))^2))), numeric(1))

dd$strategy <- with(dd, ifelse(on_front, "Best trade-off",
                               ifelse(balance > 0.65, "Middling on all axes",
                                      paste0("Held back by ", limiting_objective))))
dd$metric_set <- "von_mises_stress"

cat("\nStrategy:\n");           print(table(dd$strategy))
cat("\nLimiting objective:\n"); print(table(dd$limiting_objective))

# ── 5. SENSITIVITY TO THE STRESS OBJECTIVE ────────────────────────────────────
# Check how much the front changes if stress is left out.
stress_var <- "von_mises_stress"
nd_no_stress <- is_nd(M[, setdiff(colnames(M), stress_var), drop = FALSE])
cat(sprintf("\nFront with %s: %d | without: %d | agreement: %.1f%%\n",
            stress_var, sum(dd$on_front), sum(nd_no_stress),
            100 * mean(dd$on_front == nd_no_stress)))

# ── 6. EXPORT ─────────────────────────────────────────────────────────────────
# Save the summary table and copy the new columns back into the main data frame.
out_cols <- intersect(c("species","clade","Order","Diet.1","Diet.2","Midpoint",
                        "Period.Name","Wingspan_cm", names(OBJECTIVES),
                        "goldberg_rank","pareto_rank_ratio","on_front",
                        "dominance_count","limiting_objective",
                        "limiting_percentile","balance","distance_to_front",
                        "strategy","metric_set",
                        paste0("pareto_score_", names(OBJECTIVES)),
                        paste0("pareto_tradeoff_", TRADEOFF_CORE)), colnames(dd))

out_file <- "output/results/supplementals/wing_classification_von_mises_stress.csv"
write.csv(dd[order(-dd$pareto_rank_ratio), out_cols], out_file, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", out_file))

for (nm in c("goldberg_rank","pareto_rank_ratio","on_front","strategy",
             "limiting_objective","distance_to_front","dominance_count",
             paste0("pareto_score_", names(OBJECTIVES)),
             paste0("pareto_tradeoff_", TRADEOFF_CORE))) {
  performance_data_clean[[nm]] <- NA
  performance_data_clean[[nm]][keep] <- dd[[nm]]
}
cat(sprintf("Columns written back for %d of %d specimens\n\n", sum(keep), length(keep)))
