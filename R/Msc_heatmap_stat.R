if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")
for (dd in c("output/plots", "output/plots_PDF")) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

d <- performance_data_clean
if ("Period.Name" %in% colnames(d)) {
  pc <- trimws(as.character(d$Period.Name))
  pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
  d$Period_binned <- factor(pc, levels = c("Late Triassic", "Early+Middle Jurassic",
                                           "Late Jurassic", "Early Cretaceous",
                                           "Late Cretaceous"))
}

PAIRS <- list(
  c("Environment", "Depositional.settings.paleoenvironment"),
  c("Diet.1", "Environment"),
  c("Diet_combined", "Environment"),
  c("clade", "Diet.1"),
  c("clade", "Environment"),
  c("clade", "Period_binned"),
  c("Diet.1", "Period_binned"),
  c("Depositional.settings.paleoenvironment", "Period_binned"))

PAIRS <- Filter(function(p) all(p %in% colnames(d)), PAIRS)
cat(sprintf("Factor pairs to test: %d\n\n", length(PAIRS)))

# --------------------------------------------------------------------------------
# Cramer's V — effect size for a contingency table
# --------------------------------------------------------------------------------
# A p-value alone says whether an association exists; V says how strong it is,
# on 0 to 1, and is comparable between tables of different sizes.
cramers_v <- function(tb) {
  chi <- suppressWarnings(chisq.test(tb)$statistic)
  n <- sum(tb)
  sqrt(as.numeric(chi) / (n * (min(dim(tb)) - 1)))
}

test_rows <- list()

for (pr in PAIRS) {
  
  x <- trimws(as.character(d[[pr[1]]]))
  y <- trimws(as.character(d[[pr[2]]]))
  ok <- !is.na(x) & x != "" & !is.na(y) & y != ""
  if (sum(ok) < 20) { cat("SKIPPED", paste(pr, collapse = " x "), "- too few rows\n"); next }
  
  tb <- table(x[ok], y[ok])
  # Drop levels too rare to contribute: they add empty rows and destabilise chi2
  tb <- tb[rowSums(tb) >= MIN_N, colSums(tb) >= MIN_N, drop = FALSE]
  if (nrow(tb) < 2 || ncol(tb) < 2) {
    cat("SKIPPED", paste(pr, collapse = " x "), "- fewer than two usable levels\n"); next }
  
  ft <- fisher.test(tb, simulate.p.value = TRUE, B = 10000)
  vv <- cramers_v(tb)
  
  pair_name <- paste(gsub("\\.", "_", pr), collapse = "_x_")
  cat(sprintf("%-52s V = %.3f  p = %.4g\n", paste(pr, collapse = " x "), vv, ft$p.value))
  
  test_rows[[pair_name]] <- data.frame(
    factor_1 = pr[1], factor_2 = pr[2], n = sum(tb),
    levels_1 = nrow(tb), levels_2 = ncol(tb),
    cramers_V = round(vv, 4), fisher_p = signif(ft$p.value, 4),
    stringsAsFactors = FALSE)
  
  # --- Counts and standardised residuals -------------------------------------
  res <- suppressWarnings(chisq.test(tb)$stdres)
  
  long <- as.data.frame(as.table(tb))
  names(long) <- c("level_1", "level_2", "n")
  long$std_residual <- round(as.vector(res), 3)
  long$over_represented <- long$std_residual > 2   # |z| > 2 is the usual flag
  
  write_supp(long, paste0("S37_cooccurrence_", pair_name))
  
  # --- Figure ----------------------------------------------------------------
  p <- ggplot(long, aes(x = level_1, y = level_2, fill = std_residual)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = n), size = 2.9, colour = "grey15") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, name = "Std.\nresidual") +
    labs(title = paste(pr[1], "vs", pr[2]),
         subtitle = sprintf(
           "Cell labels: species counts. Fill: standardised residual (red = more than expected). Cramer's V = %.3f, Fisher p = %.3g",
           vv, ft$p.value),
         x = pr[1], y = pr[2]) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.subtitle = element_text(size = 7.5, colour = "grey35"),
          panel.grid = element_blank())
  
  w <- max(7, min(16, 3 + ncol(tb) * 0.9))
  h <- max(5, min(14, 3 + nrow(tb) * 0.45))
  ggsave(sprintf("output/plots/HEATMAP_%s.png", pair_name), p,width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(sprintf("output/plots_PDF/HEATMAP_%s.pdf", pair_name), p, width = w, height = h, limitsize = FALSE)
}

# --------------------------------------------------------------------------------
# Summary of all associations
# --------------------------------------------------------------------------------
if (length(test_rows) > 0) {
  tests <- do.call(rbind, test_rows)
  tests$p_adjusted_BH <- signif(p.adjust(tests$fisher_p, "BH"), 4)
  tests$strength <- cut(tests$cramers_V, c(-Inf, 0.1, 0.3, 0.5, Inf),
                        labels = c("negligible", "weak", "moderate", "strong"))
  tests <- tests[order(-tests$cramers_V), ]
  
  cat("\nASSOCIATION STRENGTH BETWEEN FACTORS:\n")
  print(tests[, c("factor_1", "factor_2", "n", "cramers_V", "p_adjusted_BH",
                  "strength")], row.names = FALSE)
  write_supp(tests, "S38_cooccurrence_tests")
  
  cat("\nStrong associations mean the two factors are confounded: any effect\n")
  cat("attributed to one may belong to the other. Watch in particular for\n")
  cat("clade x diet and setting x period - if those are strong, an ecological\n")
  cat("result cannot be separated from phylogeny or from time without a model\n")
  cat("containing both.\n\n")
}