cat("\n================================================================================\n")
cat("PHYLOMORPHOSPACE FIGURES\n")
cat("================================================================================\n\n")

FIG_W <- 9; FIG_H <- 9
# ── 4. PANEL BUILDER ──────────────────────────────────────────────────────────
ax1 <- sprintf("PC1 (%.1f%%)", summary(pca_performance)$importance[2, 1] * 100)
ax2 <- sprintf("PC2 (%.1f%%)", summary(pca_performance)$importance[2, 2] * 100)

# panel_group is rebuilt from the clade column rather than trusted from
# Msc_phylo_prep_v3.R: if its labels differ by so much as a space, the subgroup
# filter returns zero rows and the panel is drawn with branches but no tips.
if (!"clade" %in% names(pd)) stop("pd has no clade column.")
pd <- as.data.frame(pd)   # a tbl_df with list columns breaks geom_point
pd$panel_group <- ifelse(pd$clade %in% PTERODACT,
                         "Pterodactyliform", "Non-pterodactyliform")

cat("\nClades found in pd:\n"); print(table(pd$clade, useNA = "ifany"))
cat("\nPanel groups:\n");       print(table(pd$panel_group))

unknown <- setdiff(unique(pd$clade),
                   c(PTERODACT, "Darwinoptera", "Anurognathidae",
                     "Rhamphorhynchidae", "Non-breviquartossan", "Basal pterosaur"))
if (length(unknown) > 0) {
  cat("\nClades not in PTERODACT and assumed basal:\n"); print(unknown)
}

make_panel <- function(group_name, letter, subtitle) {
  
  grp <- pd[pd$panel_group == group_name, ]
  cat(sprintf("\n  %s: %d tips\n", group_name, nrow(grp)))
  if (nrow(grp) == 0) {
    cat("  WARNING: no tips in this group - landscape and grey tips only.\n")
  }
  
  p <- ggplot() +
    geom_raster(data = grid[!is.na(grid$opt), ],
                aes(x = PC1, y = PC2, fill = opt), interpolate = TRUE) +
    scale_fill_gradientn(
      colours = c("#ffffff","#e8f4ea","#bfe0c6","#8fc79c","#54a468","#1b7837"),
      name = "Pareto\noptimality", limits = c(0, 1), na.value = "white")
  
  if (!is.null(wing_perf)) {
    p <- p + geom_polygon(data = wing_perf, aes(x = PC1, y = PC2, group = grid_id),
                          fill = "grey55", colour = "grey45",
                          alpha = 0.15, linewidth = 0.18)
  }
  
  p <- p +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey35", alpha = 0.45, linewidth = 0.28) +
    # Every tip in pale grey first, so no branch ends in mid-air: edges covers
    # the whole tree while grp is only the subgroup
    geom_point(data = pd, aes(x = PC1, y = PC2),
               shape = 21, fill = "grey82", colour = "grey50",
               size = 1.6, stroke = 0.3)
  
  if (nrow(grp) > 0) {
    p <- p + geom_point(data = grp, aes(x = PC1, y = PC2),
                        shape = 21, fill = "black", colour = "black",
                        size = 1.7 , stroke = 0.4)
  }
  
  # geom_raster sets the panel limits from the grid alone; tips falling outside
  # it would be clipped and vanish without warning
  p +
    coord_fixed(ratio = 0.75, expand = FALSE,
xlim = range(c(grid$PC1, pd$PC1), na.rm = TRUE),
 ylim = range(c(grid$PC2, pd$PC2), na.rm = TRUE)) +
labs(title = sprintf("(%s) %s Phylomorphospace", letter, group_name),
         subtitle = sprintf("%s | %d of %d tips highlighted",
                            subtitle, nrow(grp), nrow(pd)),
         x = ax1, y = ax2) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 9.5, colour = "grey35"),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.6),
          panel.grid = element_blank())
}

# ── 4b. BUILD AND WRITE THE PANELS ────────────────────────────────────────────
# OUTSIDE make_panel(). Placing this loop inside the function body made it call
# itself: infinite recursion, and no figure ever written.
for (pn in list(
  list(g = "Non-pterodactyliform", l = "a",
       s = "Basal pterosaurs on the performance optimality landscape",
       f = "PERF_15A_phylo_nonpterodact"),
  list(g = "Pterodactyliform", l = "b",
       s = "Derived pterosaurs on the performance optimality landscape",
       f = "PERF_15B_phylo_pterodact"))) {
  
  p <- make_panel(pn$g, pn$l, pn$s)
  ggsave(sprintf("output/plots/%s.png", pn$f), p,
         width = FIG_W, height = FIG_H, dpi = FIG_DPI)
  ggsave(sprintf("output/plots_PDF/%s.pdf", pn$f), p,
         width = FIG_W, height = FIG_H)
  cat("Written:", pn$f, "\n")
}

# ── 5. OPTIMALITY BY GRADE - VIOLIN ───────────────────────────────────────────
vdf <- as.data.frame(pd)
vdf$grade <- ifelse(vdf$clade %in% c(PTERODACT, "Darwinoptera"),
                    "Pterodactyliformes", "Non-Pterodactyliformes")

# NA in grade or opt does not merely drop rows: a violin whose group is NA is
# drawn with the na.value colour, which is grey
vdf <- vdf[!is.na(vdf$grade) & is.finite(vdf$opt), ]

# Explicit factor levels, so ggplot cannot re-sort them out of step with the
# colour vector
vdf$grade <- factor(vdf$grade,
                    levels = c("Non-Pterodactyliformes", "Pterodactyliformes"))

cat("\nGrade levels:\n"); print(table(vdf$grade, useNA = "always"))

if (nlevels(droplevels(vdf$grade)) == 2) {
  
  GRADE_COLS <- c("Non-Pterodactyliformes" = "#4292c6",
                  "Pterodactyliformes"     = "#ef6548")
  
  wx <- wilcox.test(opt ~ grade, data = vdf)
  cat(sprintf("Wilcoxon: W = %.0f, p = %.3e\n", wx$statistic, wx$p.value))
  
  yr_v <- range(vdf$opt, na.rm = TRUE)
  
  p_violin <- ggplot(vdf, aes(x = grade, y = opt)) +
    # fill inside the geom's own aes(), not inherited from ggplot(): an
    # inherited fill is overridden by the boxplot layer, which sets fill="white"
    geom_violin(aes(fill = grade), alpha = 0.6, trim = FALSE,
                colour = "grey25", linewidth = 0.4, scale = "width") +
    geom_boxplot(width = 0.12, fill = "white", colour = "grey25",
                 alpha = 0.9, outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1.4, colour = "grey30") +
    scale_fill_manual(values = GRADE_COLS, drop = FALSE, guide = "none") +
    annotate("text", x = 1.5, y = yr_v[2] - 0.03 * diff(yr_v),
             label = sprintf("Wilcoxon\nW = %.0f\np = %.2e",
                             wx$statistic, wx$p.value),
             size = 3, hjust = 0.5, family = "mono") +
    labs(title = "(c) Pareto optimality by grade",
         subtitle = sprintf("Performance morphospace | n = %d vs %d",
                            sum(vdf$grade == levels(vdf$grade)[1]),
                            sum(vdf$grade == levels(vdf$grade)[2])),
         x = NULL, y = "Pareto rank ratio") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey35"))
  
  ggsave("output/plots/PERF_15C_optimality_violin.png",
         p_violin, width = 7, height = 8, dpi = FIG_DPI)
  ggsave("output/plots_PDF/PERF_15C_optimality_violin.pdf",
         p_violin, width = 7, height = 8)
  cat("Written: PERF_15C_optimality_violin\n")
  
} else {
  cat("Only", nlevels(droplevels(vdf$grade)), "grade(s) - violin skipped.\n")
  print(unique(vdf$clade))
}
cat("\n")