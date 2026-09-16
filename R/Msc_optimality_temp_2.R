# Add bootstrap summaries for optimality and morphospace spread.
if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")
if (!"pareto_rank_ratio" %in% colnames(performance_data_clean)) {
  stop("pareto_rank_ratio not found - run Msc_pareto_front_v3.R first.")
}

N_BOOT <- 999
set.seed(42)

d <- performance_data_clean

if (!"Time_Bin" %in% colnames(d) && "Period.Name" %in% colnames(d)) {
  pc <- trimws(as.character(d$Period.Name))
  pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
  d$Time_Bin <- factor(pc, levels = c("Late Triassic", "Early+Middle Jurassic",
                                      "Late Jurassic", "Early Cretaceous",
                                      "Late Cretaceous"))
}

################################################################################---------------------------------------------------------------------------------
# 1. Bootstrap the main scores within each group
################################################################################--
boot_mean <- function(x, n_boot = N_BOOT) {
  x <- x[is.finite(x)]
  if (length(x) < MIN_N) return(NULL)
  bs <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  ci <- quantile(bs, c(0.025, 0.975))
  c(observed = mean(x), lo = unname(ci[1]), hi = unname(ci[2]),
    boot_sd = sd(bs), n = length(x))
}

# The raw metrics replace the deleted per-objective Pareto sub-scores: they are
# what those sub-scores were derived from, on an interpretable scale.
TRACK <- intersect(c("pareto_rank_ratio", "aspect_ratio", "r2_hat", "von_mises_stress",
                     "wing_curvature", "shape_complexity", "wing_loading"),
                   colnames(d))

boot_by_group <- function(group_col) {
  if (!group_col %in% colnames(d)) return(NULL)
  g <- as.character(d[[group_col]])
  ok <- !is.na(g) & g != ""

  do.call(rbind, lapply(unique(g[ok]), function(lv) {
    s <- d[ok & g == lv, ]
    if (nrow(s) < MIN_N) return(NULL)
    do.call(rbind, lapply(TRACK, function(m) {
      r <- boot_mean(s[[m]])
      if (is.null(r)) return(NULL)
      data.frame(group = lv, metric = m, n = as.integer(r["n"]),
                 observed = round(r["observed"], 5),
                 boot_ci_lower = round(r["lo"], 5),
                 boot_ci_upper = round(r["hi"], 5),
                 boot_sd = round(r["boot_sd"], 5), stringsAsFactors = FALSE)
    }))
  }))
}

GROUPINGS <- list(
  list(col = "clade",         file = "S33_boot_optimality_by_clade"),
  list(col = "Order",         file = "S34_boot_optimality_by_order"),
  list(col = "Diet_combined", file = "S35_boot_optimality_by_diet"))

for (gg in GROUPINGS) {
  tab <- boot_by_group(gg$col)
  if (is.null(tab)) { cat("SKIPPED", gg$col, "- absent or no group reaches n >=",
                          MIN_N, "\n"); next }
  names(tab)[1] <- gg$col

  # Print the optimality rows only; the full table goes to the CSV
  opt <- tab[tab$metric == "pareto_rank_ratio", ]
  opt <- opt[order(-opt$observed), ]
  cat(sprintf("\nPareto rank ratio by %s (%d groups):\n", gg$col, nrow(opt)))
  print(opt[, c(gg$col, "n", "observed", "boot_ci_lower", "boot_ci_upper")],
        row.names = FALSE)
  write_supp(tab, gg$file)
}

# --------------------------------------------------------------------------------
# 2. Measure morphospace spread within each time bin
# --------------------------------------------------------------------------------
# Kept from the old 9B.4, which had already been corrected: distances are
# computed WITHIN each bin. Averaging per-specimen precomputed values would
# include between-bin distances and inflate the older, sparser bins.
if (all(c("PC1", "PC2") %in% colnames(d)) && "Time_Bin" %in% colnames(d)) {
  mpd_tab <- do.call(rbind, lapply(levels(d$Time_Bin), function(tb) {
    s <- d[!is.na(d$Time_Bin) & d$Time_Bin == tb &
             is.finite(d$PC1) & is.finite(d$PC2), ]
    if (nrow(s) < MIN_N) return(NULL)

    xy <- as.matrix(s[, c("PC1", "PC2")])
    obs <- mean(dist(xy))
    bs <- replicate(N_BOOT, {
      idx <- sample(nrow(xy), nrow(xy), replace = TRUE)
      mean(dist(xy[idx, , drop = FALSE]))
    })
    ci <- quantile(bs, c(0.025, 0.975))
    data.frame(Time_Bin = tb, n = nrow(s), observed_mpd = round(obs, 5),
    boot_ci_lower = round(unname(ci[1]), 5),boot_ci_upper = round(unname(ci[2]), 5),
               boot_sd = round(sd(bs), 5), stringsAsFactors = FALSE)
  }))

  if (!is.null(mpd_tab)) {
    mpd_tab$Time_Bin <- factor(mpd_tab$Time_Bin, levels = levels(d$Time_Bin))
    mpd_tab <- mpd_tab[order(mpd_tab$Time_Bin), ]
    cat("\nTrue mean pairwise distance by time bin (performance PCs):\n")
    print(mpd_tab, row.names = FALSE)
    write_supp(mpd_tab, "S36_boot_true_MPD_by_time_bin")
  }
}

cat("\nOptimality through time, with its PGLS test, is in S25/S26 from\n")
cat("Msc_temp_optimality_v1.R and is deliberately not duplicated here.\n\n")
