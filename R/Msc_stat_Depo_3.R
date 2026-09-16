cat("\n--- S39-S42: TAPHONOMIC SENSITIVITY OF THE RESULTS ---\n\n")

# Check whether preservation quality could be driving the main patterns.
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_06a_postcranial.R first.")
if (!"pareto_rank_ratio" %in% colnames(performance_data_clean)) {
  stop("pareto_rank_ratio not found - run Msc_pareto_front_v3.R first.")
}

d <- POSTCRANIAL
d$completeness <- d$n_elements / length(GM_ELEMENTS)
d <- d[is.finite(d$completeness), ]

# The threshold is arbitrary, so it is stated and its effect on sample size is
# reported. 0.8 keeps specimens with at least four fifths of elements measured.
COMPLETENESS_CUT <- 0.8
d$well_preserved <- d$completeness >= COMPLETENESS_CUT

cat(sprintf("Completeness threshold: %.2f\n", COMPLETENESS_CUT))
cat(sprintf("Well preserved: %d | Poorly preserved: %d\n\n",
            sum(d$well_preserved), sum(!d$well_preserved)))

if (sum(d$well_preserved) < 30) {
  cat("WARNING: fewer than 30 well-preserved specimens. Lower COMPLETENESS_CUT\n")
  cat("or treat the subset comparisons below as indicative only.\n\n")
}

# --------------------------------------------------------------------------------
# 1. Does preservation predict front membership?
# --------------------------------------------------------------------------------
# If it does, the "best trade-off" wings are partly the best-preserved ones, and
# the Pareto result carries a taphonomic component.
front_rows <- NULL
if ("on_front" %in% colnames(d)) {
  ok <- is.finite(d$completeness) & !is.na(d$on_front)
  tb <- table(d$well_preserved[ok], d$on_front[ok])

  if (all(dim(tb) == c(2, 2))) {
    ft <- fisher.test(tb)
    glmfit <- glm(on_front ~ completeness, data = d[ok, ], family = binomial)
    cf <- summary(glmfit)$coefficients

    front_rows <- data.frame(
      test = c("Fisher: well-preserved vs on-front",
               "Logistic: on_front ~ completeness"),
      statistic = c(round(ft$estimate, 4), round(cf[2, 1], 4)),
      p_value = signif(c(ft$p.value, cf[2, 4]), 4),
      stringsAsFactors = FALSE)

    cat("Preservation and Pareto front membership:\n")
    print(tb)
    cat(sprintf("Fisher odds ratio = %.3f, p = %.4g\n", ft$estimate, ft$p.value))
    cat(sprintf("Logistic slope on completeness = %+.4f, p = %.4g\n\n",
                cf[2, 1], cf[2, 4]))

    write_supp(front_rows, "S39_preservation_predicts_front")
  }
}

# --------------------------------------------------------------------------------
# 2. Do the group tests give the same verdict on the well-preserved subset?
# --------------------------------------------------------------------------------
# The core sensitivity analysis: every Kruskal-Wallis is run twice, once on the
# full sample and once on the well-preserved subset, and the verdicts compared.
TEST_METRICS <- intersect(METRICS, colnames(d))
TEST_GROUPS  <- intersect(c("clade", "Order", "Diet.1", "Diet_combined"),
                          colnames(d))

kw_p <- function(sub, m, g) {
  s <- sub[is.finite(sub[[m]]) & !is.na(sub[[g]]) & sub[[g]] != "", ]
  keep_lv <- names(which(table(s[[g]]) >= MIN_N))
  s <- s[s[[g]] %in% keep_lv, ]
  if (length(keep_lv) < 2 || nrow(s) < 15) return(NA_real_)
  kruskal.test(as.formula(paste0("`", m, "` ~ `", g, "`")), data = s)$p.value
}

sens <- do.call(rbind, lapply(TEST_GROUPS, function(g) {
  do.call(rbind, lapply(TEST_METRICS, function(m) {
    p_all  <- kw_p(d, m, g)
    p_well <- kw_p(d[d$well_preserved, ], m, g)
    if (is.na(p_all) && is.na(p_well)) return(NULL)
    data.frame(group = g, metric = m,
               n_all = sum(is.finite(d[[m]])),
               n_well = sum(is.finite(d[[m]]) & d$well_preserved),
               p_full_sample = signif(p_all, 4),
               p_well_preserved = signif(p_well, 4),
               stringsAsFactors = FALSE)
  }))
}))

if (!is.null(sens)) {
  sens$verdict <- with(sens, ifelse(
    is.na(p_full_sample) | is.na(p_well_preserved), "not testable",
    ifelse(p_full_sample < 0.05 & p_well_preserved < 0.05, "robust",
           ifelse(p_full_sample >= 0.05 & p_well_preserved >= 0.05, "null in both",
                  ifelse(p_full_sample < 0.05, "LOST in well-preserved subset",
                         "appears only in well-preserved subset")))))

  cat("Sensitivity of the group tests to preservation:\n")
  print(table(sens$verdict))
  fragile <- sens[grepl("LOST|appears only", sens$verdict), ]
  if (nrow(fragile) > 0) {
    cat("\nResults that change with the subset:\n")
    print(fragile[, c("group", "metric", "p_full_sample", "p_well_preserved",
                      "verdict")], row.names = FALSE)
  }
  write_supp(sens, "S40_subset_sensitivity")
}

# --------------------------------------------------------------------------------
# 3. Partial correlations, controlling for completeness
# --------------------------------------------------------------------------------
# This asks whether a correlation stays once preservation quality is held constant.
# A metric-size relationship that vanishes once completeness is held constant
# was a preservation artefact. Computed as the correlation of the residuals of
# both variables on completeness.
partial_cor <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 20) return(c(NA, NA, NA))
  rx <- residuals(lm(x[ok] ~ z[ok]))
  ry <- residuals(lm(y[ok] ~ z[ok]))
  raw <- cor(x[ok], y[ok], method = "spearman")
  ct <- cor.test(rx, ry, method = "spearman", exact = FALSE)
  c(raw, unname(ct$estimate), ct$p.value)
}

SIZE_REF <- intersect(c("Wingspan_cm", "SkeletalSize"), colnames(d))[1]

pc_tab <- NULL
if (!is.na(SIZE_REF)) {
  pc_tab <- do.call(rbind, lapply(TEST_METRICS, function(m) {
    r <- partial_cor(d[[m]], d[[SIZE_REF]], d$completeness)
    if (all(is.na(r))) return(NULL)
    data.frame(metric = m, reference = SIZE_REF,
               rho_raw = round(r[1], 4),
               rho_partial = round(r[2], 4),
               p_partial = signif(r[3], 4),
               change = round(r[2] - r[1], 4), stringsAsFactors = FALSE)
  }))

  if (!is.null(pc_tab)) {
    pc_tab$artefact_suspected <- abs(pc_tab$change) > 0.1
    cat("\nCorrelation with", SIZE_REF, "before and after controlling for completeness:\n")
    print(pc_tab, row.names = FALSE)
    write_supp(pc_tab, "S41_partial_correlations")
  }
}
