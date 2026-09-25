# Test whether groups occupy different parts of the multivariate performance space.
if (!exists("MV_MATRIX")) stop("Source Msc_stats_00_setup.R first.")
if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required: install.packages(\"vegan\")")
}
library(vegan)

cat(sprintf("Matrix: %d specimens x %d metrics | %d permutations\n",
            nrow(MV_MATRIX), ncol(MV_MATRIX), N_PERM))
cat("Factors tested:", paste(FACTORS, collapse = ", "), "\n")
if ("Diet_combined" %in% FACTORS) {
  n_lv <- sum(table(MV_META$Diet_combined) >= MIN_N)
  cat(sprintf("Diet_combined: %d levels reach n >= %d\n", n_lv, MIN_N))
}
cat("\n")

# Turn the multivariate table into distances between specimens.
D_FULL <- vegdist(MV_MATRIX, method = "euclidean")

perm_rows <- list(); disp_rows <- list()

for (f in FACTORS) {
  grp <- as.character(MV_META[[f]])
  ok <- !is.na(grp) & grp != ""
  # Drop groups that are too small for a meaningful permutation test.
  keep_lv <- names(which(table(grp[ok]) >= MIN_N))
  ok <- ok & grp %in% keep_lv

  if (sum(ok) < 20 || length(unique(grp[ok])) < 2) {
    cat("SKIPPED", f, "- too few usable specimens or levels\n"); next
  }

  Ds <- as.dist(as.matrix(D_FULL)[ok, ok])
  g  <- factor(grp[ok])

  ad <- adonis2(Ds ~ g, permutations = N_PERM)
  bd <- betadisper(Ds, g)
  bp <- permutest(bd, permutations = N_PERM)

  perm_rows[[f]] <- data.frame(
    factor = f, n = sum(ok), levels = nlevels(g), df = ad$Df[1],
    SumOfSqs = round(ad$SumOfSqs[1], 3), R2 = round(ad$R2[1], 4),
    F_value = round(ad$F[1], 3), p_value = ad$`Pr(>F)`[1],
    dispersion_F = round(bp$tab$F[1], 3), dispersion_p = bp$tab$`Pr(>F)`[1],
    stringsAsFactors = FALSE)

  disp_rows[[f]] <- data.frame(
    factor = f, group = levels(g), n = as.integer(table(g)),
    mean_distance_to_centroid = round(as.numeric(tapply(bd$distances, g, mean)), 4),
    sd_distance = round(as.numeric(tapply(bd$distances, g, sd)), 4),
    stringsAsFactors = FALSE)

  cat(sprintf("%-42s R2=%.3f  p=%.4g  | dispersion p=%.4g\n",
              f, ad$R2[1], ad$`Pr(>F)`[1], bp$tab$`Pr(>F)`[1]))
}

if (length(perm_rows) == 0) stop("No factor produced a usable PERMANOVA.")

perm_tab <- do.call(rbind, perm_rows)
perm_tab$p_adjusted_BH <- signif(p.adjust(perm_tab$p_value, method = "BH"), 4)

# Add a plain-language verdict to each row so the main result is easy to read.
# This matters because a spread effect is not the same as a location effect.
perm_tab$interpretation <- ifelse(
  perm_tab$p_adjusted_BH >= 0.05, "no group effect",
  ifelse(perm_tab$dispersion_p < 0.05,
         "groups differ, but dispersion differs too - not a pure location effect",
         "groups differ in location"))

cat("\nPERMANOVA SUMMARY:\n")
# Show a fixed number of decimals so every row lines up.
perm_tab_display <- perm_tab
perm_tab_display$R2             <- format_fixed(perm_tab$R2, 4)
perm_tab_display$F_value        <- format_fixed(perm_tab$F_value, 3)
perm_tab_display$p_adjusted_BH  <- format_fixed(perm_tab$p_adjusted_BH, 4)
perm_tab_display$dispersion_p   <- format_fixed(perm_tab$dispersion_p, 4)
print(perm_tab_display[, c("factor", "n", "levels", "R2", "F_value",
                   "p_adjusted_BH", "dispersion_p", "interpretation")],
      row.names = FALSE)

disp_tab <- do.call(rbind, disp_rows)
disp_tab$mean_distance_to_centroid <- format_fixed(disp_tab$mean_distance_to_centroid, 4)
disp_tab$sd_distance                <- format_fixed(disp_tab$sd_distance, 4)

write_supp(perm_tab_display, "S6_permanova_summary")
write_supp(disp_tab, "S7_dispersion_by_group")

cat("\nRead every clade result against the phylANOVA in Msc_phylo_signal_v3.R:\n")
cat("clades are defined on the tree, so a significant clade effect is expected\n")
cat("regardless of function.\n\n")
