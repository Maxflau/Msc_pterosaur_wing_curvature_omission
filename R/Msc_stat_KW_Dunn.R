if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")

# Run the group comparison tests and save their output tables.
# --------------------------------------------------------------------------------
# 1. Kruskal-Wallis, one test per metric per taxonomic level
# --------------------------------------------------------------------------------
TEST_LEVELS <- intersect(c(GROUPS, "Diet_combined"),
                         colnames(performance_data_clean))

for (g in TEST_LEVELS) {
  d <- performance_data_clean[!is.na(performance_data_clean[[g]]), ]

  kw <- do.call(rbind, lapply(METRICS, function(m) {
    sub <- d[is.finite(d[[m]]), c(g, m)]
    if (length(unique(sub[[g]])) < 2 || length(unique(sub[[m]])) < 2) return(NULL)
    k <- kruskal.test(as.formula(paste0("`", m, "` ~ `", g, "`")), data = sub)
    data.frame(metric = m,
               chi_squared = round(unname(k$statistic), 4),
               df = unname(k$parameter),
               n = nrow(sub),
               p_value = signif(k$p.value, 4),
               stringsAsFactors = FALSE)
  }))

  if (is.null(kw)) { cat("SKIPPED", g, "- no testable metric\n"); next }

  kw$p_adjusted_BH <- signif(p.adjust(kw$p_value, method = "BH"), 4)
  kw$significant   <- kw$p_adjusted_BH < 0.05

  cat(sprintf("Kruskal-Wallis by %s (%d metrics, %d significant after BH):\n",
              g, nrow(kw), sum(kw$significant, na.rm = TRUE)))
  # Show a fixed number of decimals so every row lines up (e.g. "0.0500", not ".05").
  kw_display <- kw
  kw_display$chi_squared   <- format_fixed(kw$chi_squared, 4)
  kw_display$p_value       <- format_fixed(kw$p_value, 4)
  kw_display$p_adjusted_BH <- format_fixed(kw$p_adjusted_BH, 4)
  print(kw_display, row.names = FALSE)
  write_supp(kw_display, paste0("S2_kruskal_wallis_by_", tolower(g)))
  cat("\n")
}

# --------------------------------------------------------------------------------
# 2. Dunn pairwise comparisons
# --------------------------------------------------------------------------------
# If a global test is possible, compare the usable groups two at a time.
if (!requireNamespace("FSA", quietly = TRUE)) {
  cat("Package 'FSA' not installed - pairwise Dunn tests skipped.\n")
  cat("Run: install.packages(\"FSA\")\n\n")
} else {
  for (g in TEST_LEVELS) {
    d <- performance_data_clean[!is.na(performance_data_clean[[g]]), ]
    keep_lv <- names(which(table(d[[g]]) >= MIN_N))
    d <- d[d[[g]] %in% keep_lv, ]

    if (length(keep_lv) < 2) { cat("SKIPPED", g, "- fewer than two usable groups\n"); next }

    pw <- do.call(rbind, lapply(METRICS, function(m) {
      sub <- d[is.finite(d[[m]]), c(g, m)]
      if (length(unique(sub[[g]])) < 2) return(NULL)
      res <- tryCatch(
        FSA::dunnTest(reformulate(g, response = m), data = sub, method = "bh")$res,
        error = function(e) NULL)
      if (is.null(res)) return(NULL)
      parts <- strsplit(as.character(res$Comparison), " - ")
      data.frame(metric = m,
                 group_1 = trimws(vapply(parts, function(p) p[1], character(1))),
                 group_2 = trimws(vapply(parts, function(p) p[2], character(1))),
                 Z = round(res$Z, 3),
                 p_unadjusted = signif(res$P.unadj, 4),
                 p_adjusted_BH = signif(res$P.adj, 4),
                 significant = res$P.adj < 0.05,
                 stringsAsFactors = FALSE)
    }))

    if (is.null(pw)) { cat("SKIPPED", g, "- no pairwise result\n"); next }

    cat(sprintf("Dunn pairwise by %s: %d groups, %d comparisons, %d significant\n",
                g, length(keep_lv), nrow(pw), sum(pw$significant, na.rm = TRUE)))
    # Show a fixed number of decimals so every row lines up.
    pw_display <- pw
    pw_display$Z             <- format_fixed(pw$Z, 3)
    pw_display$p_unadjusted  <- format_fixed(pw$p_unadjusted, 4)
    pw_display$p_adjusted_BH <- format_fixed(pw$p_adjusted_BH, 4)
    write_supp(pw_display, paste0("S3_pairwise_dunn_by_", tolower(g)))
  }
}

cat("\n")
