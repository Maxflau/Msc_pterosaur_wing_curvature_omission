cat("\n================================================================================\n")
cat("PHYLOMORPHOSPACE FIGURES\n")
cat("================================================================================\n\n")

library(ggplot2)

# Make sure the plot folders exist before saving figures.
for (dd in c("output/plots", "output/plots_PDF")) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}
if (!exists("phy_pruned") || !exists("phylo_data")) {
  stop("phy_pruned / phylo_data not found - source Msc_phylo_prep_v3.R first.")
}

# Use the Pareto score already stored in the main data table.
OPT_VAR <- if ("pareto_rank_ratio" %in% colnames(performance_data_clean)) {
  "pareto_rank_ratio"
} else stop("pareto_rank_ratio not found - run Msc_pareto_front_v3.R first.")

FIG_W <- 12; FIG_H <- 8; FIG_DPI <- 300

PTERODACT <- c("Azhdarchoidea", "Ctenochasmatoidea", "Dsungaripteroidea",
               "Ornithocheiromorpha", "Pteranodontia", "Basal pterodactyliform")

# ── 1. ANCESTRAL STATES ON THE PERFORMANCE AXES ───────────────────────────────
# Rebuild the tree coordinates so the branches match the plotted points.
pd <- phylo_data[is.finite(phylo_data$PC1) & is.finite(phylo_data$PC2), ]
tr <- drop.tip(phy_pruned,
               setdiff(phy_pruned$tip.label,
                       phy_pruned$tip.label[normalise_name(phy_pruned$tip.label)
                                            %in% normalise_name(pd$species)]))
pd <- pd[match(normalise_name(tr$tip.label), normalise_name(pd$species)), ]
cat(sprintf("Tips on the performance plane: %d\n", nrow(pd)))

s1 <- setNames(pd$PC1, tr$tip.label); s2 <- setNames(pd$PC2, tr$tip.label)
a1 <- fastAnc(tr, s1); a2 <- fastAnc(tr, s2)
all1 <- c(s1, a1); all2 <- c(s2, a2)

nt <- Ntip(tr)
child_nm <- ifelse(tr$edge[, 2] <= nt, tr$tip.label[tr$edge[, 2]],
                   as.character(tr$edge[, 2]))
edges <- data.frame(x = all1[as.character(tr$edge[, 1])],
                    y = all2[as.character(tr$edge[, 1])],
                    xend = all1[child_nm], yend = all2[child_nm])
edges <- edges[complete.cases(edges), ]

k1 <- phylosig(tr, s1, method = "K", test = TRUE, nsim = 1000)
k2 <- phylosig(tr, s2, method = "K", test = TRUE, nsim = 1000)
k_text <- sprintf("Blomberg's K\nPC1: K=%.3f  p=%.3f\nPC2: K=%.3f  p=%.3f",
                  k1$K, k1$P, k2$K, k2$P)
cat(k_text, "\n\n")

# ── 2. OPTIMALITY SURFACE ─────────────────────────────────────────────────────
# Fill a grid with smoothed Pareto values across the performance space.
pd$opt <- performance_data_clean[[OPT_VAR]][match(pd$species,
                                                  performance_data_clean$species)]

r1 <- range(pd$PC1); r2 <- range(pd$PC2)
e1 <- diff(r1) * 0.12; e2 <- diff(r2) * 0.12
grid <- expand.grid(PC1 = seq(r1[1] - e1, r1[2] + e1, length.out = 90),
                    PC2 = seq(r2[1] - e2, r2[2] + e2, length.out = 90))

# Set the smoothing width from the size of the morphospace itself.
bw <- 0.10 * max(diff(r1), diff(r2))
ok <- is.finite(pd$opt)
grid$opt <- vapply(seq_len(nrow(grid)), function(i) {
  d2 <- (pd$PC1[ok] - grid$PC1[i])^2 + (pd$PC2[ok] - grid$PC2[i])^2
  w <- exp(-d2 / (2 * bw^2))
  if (sum(w) > 1e-8) weighted.mean(pd$opt[ok], w) else NA_real_
}, numeric(1))

# ── 3. WING LATTICE ACROSS THE PANEL ──────────────────────────────────────────
# Add nearby wing outlines to show what the morphospace looks like in shape terms.
wing_perf <- NULL
if (exists("outlines_list") && length(outlines_list) > 0) {
  gx <- seq(r1[1] - e1, r1[2] + e1, length.out = 13)
  gy <- seq(r2[1] - e2, r2[2] + e2, length.out = 13)
  cells <- expand.grid(x = gx, y = gy)
  sc <- 0.025 * max(diff(range(gx)), diff(range(gy)))
  out <- list()
  for (i in seq_len(nrow(cells))) {
    d2 <- (pd$PC1 - cells$x[i])^2 + (pd$PC2 - cells$y[i])^2
    o <- outlines_list[[pd$species[which.min(d2)]]]
    if (is.null(o) || nrow(o) < 3) next
    xx <- o$x - mean(o$x); yy <- o$y - mean(o$y); m <- max(abs(c(xx, yy)))
    if (!is.finite(m) || m == 0) next
    out[[length(out) + 1]] <- data.frame(PC1 = cells$x[i] + sc * xx / m,
                                         PC2 = cells$y[i] + sc * yy / m,
                                         grid_id = i)
  }
  if (length(out) > 0) wing_perf <- do.call(rbind, out)
}
cat(sprintf("Wing lattice: %d profiles\n\n",
 if (is.null(wing_perf)) 0 else length(unique(wing_perf$grid_id))))
