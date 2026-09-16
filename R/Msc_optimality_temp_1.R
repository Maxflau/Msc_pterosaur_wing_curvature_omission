# Summarise how the optimality score changes through time.
cat("\n--- S25/S26: OPTIMALITY THROUGH TIME ---\n\n")

if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")
if (!"pareto_rank_ratio" %in% colnames(performance_data_clean)) {
  stop("pareto_rank_ratio not found - run Msc_pareto_front_v3.R first.")
}

d <- performance_data_clean[is.finite(performance_data_clean$pareto_rank_ratio), ]

# Rebuild the time bins here so this script can run on its own.
if (!"Time_Bin" %in% colnames(d) && "Period.Name" %in% colnames(d)) {
  pc <- trimws(as.character(d$Period.Name))
  pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
  d$Time_Bin <- factor(pc, levels = c("Late Triassic", "Early+Middle Jurassic",
                                      "Late Jurassic", "Early Cretaceous",
                                      "Late Cretaceous"))
}
d <- d[!is.na(d$Time_Bin), ]
cat(sprintf("Specimens with an optimality score and a time bin: %d\n", nrow(d)))

# --------------------------------------------------------------------------------
# 1. Estimate the average score in each time bin
# --------------------------------------------------------------------------------
boot_ci <- function(x, n_boot = 2000) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(NA, NA, NA))
  bm <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  c(mean(x), unname(quantile(bm, 0.025)), unname(quantile(bm, 0.975)))
}

time_tab <- do.call(rbind, lapply(levels(d$Time_Bin), function(tb) {
  s <- d[d$Time_Bin == tb, ]
  if (nrow(s) == 0) return(NULL)
  ci <- boot_ci(s$pareto_rank_ratio)
  data.frame(Time_Bin = tb, n = nrow(s),
             mean_age = round(mean(s$Midpoint, na.rm = TRUE), 1),
             mean_PRR = round(ci[1], 4), lo = round(ci[2], 4), hi = round(ci[3], 4),
             pct_on_front = round(100 * mean(s$on_front, na.rm = TRUE), 1),
             stringsAsFactors = FALSE)
}))

cat("\nOptimality by time bin:\n"); print(time_tab, row.names = FALSE)
write_supp(time_tab, "S25_optimality_through_time")

# --------------------------------------------------------------------------------
# 2. Break each time bin down by clade
# --------------------------------------------------------------------------------
# A bin-level mean can rise simply because a high-scoring clade radiates. This
# table shows whether the rise happens WITHIN clades or only between them.
if ("clade" %in% colnames(d)) {
  by_cb <- do.call(rbind, lapply(levels(d$Time_Bin), function(tb) {
    s <- d[d$Time_Bin == tb, ]
    do.call(rbind, lapply(unique(s$clade), function(cl) {
      x <- s$pareto_rank_ratio[s$clade == cl]
      if (length(x) < MIN_N) return(NULL)
      data.frame(Time_Bin = tb, clade = cl, n = length(x),
                 mean_PRR = round(mean(x), 4), stringsAsFactors = FALSE)
    }))
  }))
  cat(sprintf("\nClade-by-bin cells with n >= %d: %d\n", MIN_N, nrow(by_cb)))
  write_supp(by_cb, "S25b_optimality_by_clade_and_bin")
}

# --------------------------------------------------------------------------------
# 3. Compare the plain trend and the phylogenetic trend
# --------------------------------------------------------------------------------
# Reported as a pair so the write-up cannot quote one without the other.
if ("Midpoint" %in% colnames(d)) {
  ols <- lm(pareto_rank_ratio ~ Midpoint, data = d)
  o <- summary(ols)$coefficients
  cat(sprintf("\nOLS   : slope = %+.5f, p = %.4g, R2 = %.3f\n",
              o[2, 1], o[2, 4], summary(ols)$r.squared))

  rows <- data.frame(model = "OLS", n = nrow(d), slope = round(o[2, 1], 5),
                     se = round(o[2, 2], 5), p_value = signif(o[2, 4], 4),
                     lambda = NA, stringsAsFactors = FALSE)

  if (exists("phy_pruned") && requireNamespace("nlme", quietly = TRUE) &&
      exists("normalise_name")) {
    dd <- d[normalise_name(d$species) %in% normalise_name(phy_pruned$tip.label), ]
    tr <- ape::drop.tip(phy_pruned,
                        phy_pruned$tip.label[!normalise_name(phy_pruned$tip.label)
                                             %in% normalise_name(dd$species)])
    dd <- dd[match(normalise_name(tr$tip.label), normalise_name(dd$species)), ]
    dd <- dd[complete.cases(dd[, c("pareto_rank_ratio", "Midpoint")]), ]

    m <- tryCatch(nlme::gls(pareto_rank_ratio ~ Midpoint, data = dd,
                            correlation = ape::corPagel(1, phy = tr, form = ~1),
                            method = "ML"),
                  error = function(e) { cat("PGLS failed:", e$message, "\n"); NULL })
    if (!is.null(m)) {
      tt <- summary(m)$tTable
      lam <- as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))
      cat(sprintf("PGLS  : slope = %+.5f, p = %.4g, lambda = %.3f\n",
                  tt["Midpoint", "Value"], tt["Midpoint", "p-value"], lam))
      rows <- rbind(rows, data.frame(
        model = "PGLS (Pagel lambda)", n = nrow(dd),
        slope = round(tt["Midpoint", "Value"], 5),
        se = round(tt["Midpoint", "Std.Error"], 5),
        p_value = signif(tt["Midpoint", "p-value"], 4),
        lambda = round(lam, 3), stringsAsFactors = FALSE))
    }
  } else {
    cat("PGLS skipped - phy_pruned or nlme unavailable.\n")
  }

  write_supp(rows, "S26_optimality_pgls")

  if (nrow(rows) == 2 && rows$p_value[1] < 0.05 && rows$p_value[2] >= 0.05) {
    cat("\nThe trend does NOT survive phylogenetic correction. It reflects the\n")
    cat("pterodactyloid radiation, not evolution in performance, and must not be\n")
    cat("reported as a temporal trend.\n")
  }
}
cat("\n")
