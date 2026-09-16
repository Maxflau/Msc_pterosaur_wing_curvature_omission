if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")
if (!exists("outlines_list")) stop("outlines_list not found.")

for (d in c("clade_analyses_v4/plots", "clade_analyses_v4/plots_PDF")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

#####################################################################################
# 1. Distribution of measured curvature
#####################################################################################
d <- performance_data_clean[is.finite(performance_data_clean$wing_curvature), ]
cat(sprintf("Specimens with a curvature value: %d of %d\n",
            nrow(d), nrow(performance_data_clean)))
print(summary(d$wing_curvature))

CURV_LEVELS <- intersect(c("clade", "Order", "Diet_combined"), colnames(d))

dist_tab <- do.call(rbind, lapply(CURV_LEVELS, function(g) {
  s <- d[!is.na(d[[g]]) & d[[g]] != "", ]
  do.call(rbind, lapply(unique(s[[g]]), function(lv) {
    x <- s$wing_curvature[s[[g]] == lv]
    if (length(x) < MIN_N) return(NULL)
    data.frame(level = g, group = lv, n = length(x),
               mean = round(mean(x), 5), sd = round(sd(x), 5),
               median = round(median(x), 5),
               q05 = round(quantile(x, 0.05), 5),
               q95 = round(quantile(x, 0.95), 5), stringsAsFactors = FALSE)
  }))
}))
write_supp(dist_tab, "S23_curvature_distribution")

##########################################################################################
# 2. Curvature-free baseline: the convex hull of each outline
#############################################
# chull() removes every concavity, so the hull IS a curvature-free
# reconstruction of the same wing. Measuring both with the same functions makes
# the effect of curvature a measured quantity rather than an assumption.
hull_metrics <- function(o) {
  h <- chull(o$x, o$y)
  oh <- data.frame(x = o$x[h], y = o$y[h])
  oh <- rbind(oh, oh[1, ])
  c(curv = calculate_wing_curvature(oh),
    r2   = calculate_r2_hat(oh),
    cplx = calculate_shape_complexity(oh),
    ar   = calculate_aspect_ratio(oh))
}

cat("\nMeasuring the convex-hull (curvature-free) version of each outline...\n")
hm <- t(vapply(performance_data_clean$species, function(s) {
  o <- outlines_list[[s]]
  if (is.null(o) || nrow(o) < 4) return(c(curv = NA, r2 = NA, cplx = NA, ar = NA))
  tryCatch(hull_metrics(o), error = function(e) c(curv = NA, r2 = NA, cplx = NA, ar = NA))
}, numeric(4)))

cmp <- data.frame(
  species = performance_data_clean$species,
  clade = performance_data_clean$clade,
  curvature_measured = performance_data_clean$wing_curvature,
  curvature_hull = hm[, "curv"],
  r2_measured = performance_data_clean$r2_hat, r2_hull = hm[, "r2"],
  complexity_measured = performance_data_clean$shape_complexity,
  complexity_hull = hm[, "cplx"], stringsAsFactors = FALSE)

cmp$curvature_lost <- cmp$curvature_measured - cmp$curvature_hull
cmp$r2_shift_pct <- 100 * (cmp$r2_hull - cmp$r2_measured) / cmp$r2_measured

ok <- is.finite(cmp$curvature_lost)
cat(sprintf("\nCurvature removed by hulling: median %.4f, max %.4f (n = %d)\n",
            median(cmp$curvature_lost[ok]), max(cmp$curvature_lost[ok]), sum(ok)))
cat(sprintf("Resulting shift in r2_hat: median %.2f%%, max %.2f%%\n",
            median(abs(cmp$r2_shift_pct), na.rm = TRUE),
            max(abs(cmp$r2_shift_pct), na.rm = TRUE)))

# If the outlines already come from chull() in Msc_outline_extract_1.R, the two
# columns are identical and the difference is zero everywhere. That is itself
# the finding: the current extraction discards curvature before it is measured.
if (all(abs(cmp$curvature_lost[ok]) < 1e-8)) {
  cat("\nWARNING: hull and measured curvature are identical for every specimen.\n")
  cat("The outlines are ALREADY convex hulls (Msc_outline_extract_1.R applies\n")
  cat("chull), so wing_curvature measures the bulge of a hull, not membrane\n")
  cat("camber. This must be stated in the methods.\n")
}

write_supp(cmp, "S24_curvature_free_comparison")

##########################################################################################
# 3. Figures
##########################################################################################
p1 <- ggplot(d, aes(x = wing_curvature)) +
  geom_histogram(bins = 30, fill = "#2E8B57", colour = "grey30", linewidth = 0.25) +
  geom_vline(xintercept = median(d$wing_curvature), linetype = "dashed",
             colour = "grey20") +
  labs(title = "Distribution of wing curvature across reconstructions",
       subtitle = sprintf("n = %d; dashed line: median %.4f",
                          nrow(d), median(d$wing_curvature)),
       x = "Wing curvature (camber ratio)", y = "Specimens") + theme_bw(base_size = 12) + theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))

ggsave("clade_analyses_v4/plots/CURVATURE_01_distribution.png", p1, width = 8, height = 5, dpi = 300)
ggsave("clade_analyses_v4/plots_PDF/CURVATURE_01_distribution.pdf", p1,
       width = 8, height = 5)

if ("clade" %in% colnames(d)) {
  p2 <- ggplot(d[!is.na(d$clade), ], aes(x = reorder(clade, wing_curvature,
                                                     FUN = median),
                                         y = wing_curvature, fill = clade)) +
    geom_violin(alpha = 0.5, scale = "width", colour = "grey35", linewidth = 0.3) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.5, colour = "grey20") +
    scale_fill_viridis_d(option = "mako", guide = "none", begin = 0.2, end = 0.9) +
    coord_flip() +
    labs(title = "Wing curvature by clade", x = NULL,
         y = "Wing curvature (camber ratio)") +
    theme_bw(base_size = 12) +
    theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))
  
  ggsave("clade_analyses_v4/plots/CURVATURE_02_free_vs_measured.png", p2, width = 9, height = 6, dpi = 300)
  ggsave("clade_analyses_v4/plots_PDF/CURVATURE_02_free_vs_measured.pdf", p2,width = 9, height = 6)
}