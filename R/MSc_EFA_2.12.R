# Save loadings, variance, and test tables for the EFA shape space.
bio_present <- BIO_VARS[BIO_VARS %in% names(sp)]
bio_scaled  <- scale(sp[, bio_present])

# NOTE: this is a correlation table (biomechanical metrics ~ shapePC1/PC2),
# NOT the EFA PCA rotation matrix — kept for the biplot arrows below.
efa_loadings <- data.frame(
  variable = bio_present,
  shapePC1 = apply(bio_scaled, 2, function(x) cor(x, sp$shapePC1, use="complete.obs")),
  shapePC2 = apply(bio_scaled, 2, function(x) cor(x, sp$shapePC2, use="complete.obs"))
)
save_csv(efa_loadings, "EFA_shape_pc_biomech_correlations")
cat("  EFA shapePC ~ biomechanical metric correlations:\n"); print(efa_loadings)

# ── True EFA PCA loadings (rotation matrix) ────────────────────────────────────
# Save the rotation values that define the EFA shape axes.
# Source: pca_shape object from MSc_rep.R (prcomp on z-scored EFA coefficients)
# Matches the format of the biomechanical 00_pca_loadings_recalculated.csv
# (variable, PC1, PC2) but for the EFA harmonic-coefficient PCA.
if (exists("pca_shape")) {
  rot <- unclass(pca_shape$rotation)
  pc_cols <- intersect(c("PC1","PC2"), colnames(rot))
  efa_pca_loadings <- data.frame(variable = rownames(rot), rot[, pc_cols, drop=FALSE],
                                 row.names = NULL)
  save_csv(efa_pca_loadings, "EFA_pca_loadings")
  cat("  EFA PCA loadings (rotation matrix):\n"); print(efa_pca_loadings)
} else {
  cat("  pca_shape not found — EFA_pca_loadings skipped (run MSc_rep.R first).\n")
}

# ── Variance explained by each EFA shape PC axis ──────────────────────────────
# Save how much of the total variation each axis explains.
# Source: pca_shape object from MSc_rep.R (prcomp on z-scored EFA coefficients)
# Column format matches the biomechanical 01_pca_variance_explained.csv
# (PC, Var_pct, Cumul_pct — percentages, not fractions).
if (exists("pca_shape")) {
  importance <- summary(pca_shape)$importance
  var_exp_df <- data.frame(
    PC        = colnames(importance),
    Var_pct   = round(importance["Proportion of Variance",] * 100, 2),
    Cumul_pct = round(importance["Cumulative Proportion",]  * 100, 2)
  )
  save_csv(var_exp_df, "EFA_pca_variance_explained")
  cat(sprintf("  EFA PCA: %d axes | PC1=%.1f%% | PC2=%.1f%% | Cumul PC1+PC2=%.1f%%\n",
              nrow(var_exp_df),
              var_exp_df$Var_pct[1],
              var_exp_df$Var_pct[2],
              var_exp_df$Cumul_pct[2]))
} else if (exists("shape_var")) {
  # Fallback: use shape_var vector if pca_shape not in environment
  # (shape_var is assumed to already be in percentage units)
  var_exp_df <- data.frame(
    PC        = paste0("PC", seq_along(shape_var)),
    Var_pct   = round(shape_var, 2),
    Cumul_pct = round(cumsum(shape_var), 2)
  )
  save_csv(var_exp_df, "EFA_pca_variance_explained")
  cat("  EFA variance explained saved from shape_var.\n")
} else {
  cat("  pca_shape and shape_var not found — run MSc_rep.R first.\n")
}

# Biplot
# Draw a biplot so the main directions can be seen.
arrow_s <- max(abs(sp[,c("shapePC1","shapePC2")])) * 0.55
p_biplot <- ggplot() +
  geom_point(data=sp, aes(shapePC1, shapePC2, colour=clade),
             size=2, alpha=0.75) +
  geom_segment(data=efa_loadings,
               aes(x=0,y=0, xend=shapePC1*arrow_s, yend=shapePC2*arrow_s),
               arrow=arrow(length=unit(0.2,"cm")), colour="black", linewidth=0.7) +
  geom_label_repel(data=efa_loadings,
                   aes(x=shapePC1*arrow_s, y=shapePC2*arrow_s, label=variable),
                   size=3, max.overlaps=20) +
  labs(title="EFA shape morphospace — biplot with biomechanical loadings",
       x=efa_lab1, y=efa_lab2, colour="Clade") +
  theme_bw()
save_plot(p_biplot, "EFA_biplot_loadings")

# =============================================================================
# SECTION C — PERMANOVA ON EFA SHAPE PC SPACE 
# =============================================================================
# Test how strongly the main groups differ in EFA shape space.
cat("\n── Section C: PERMANOVA on EFA shape PC space ──\n")

shape_dist <- dist(sp[,c("shapePC1","shapePC2")], method="euclidean")
set.seed(42)

tidy_perm <- function(res, factor) {
  tbl <- as.data.frame(res); tbl$Term <- rownames(tbl); tbl$Factor <- factor
  tbl <- tbl[!is.na(tbl$F),]
  names(tbl)[names(tbl)=="Pr(>F)"] <- "p_value"
  tbl %>% mutate(sig=case_when(p_value<0.001~"***",p_value<0.01~"**",
                               p_value<0.05~"*",p_value<0.1~".",TRUE~""))
}

# Helper: clean pairwise adonis2 output to a tidy data frame
tidy_pairwise <- function(pw_obj, group_label) {
  pw_rows <- lapply(names(pw_obj)[names(pw_obj)!="parent_call"], function(nm) {
    tbl <- as.data.frame(pw_obj[[nm]]); tbl$Comparison <- nm; tbl$Term <- rownames(tbl); tbl
  })
  do.call(rbind, pw_rows) %>%
    filter(!Term %in% c("Residual","Total")) %>%
    select(Comparison,Df,SumOfSqs,R2,F,`Pr(>F)`) %>%
    rename(p_adj=`Pr(>F)`) %>%
    mutate(Group=group_label,
           sig=case_when(p_adj<0.001~"***",p_adj<0.01~"**",
                         p_adj<0.05~"*",  p_adj<0.1~".",TRUE~"")) %>%
    rename(`R²`=R2)
}

avail_groups <- GROUPS[GROUPS %in% names(sp)]
# Add Time_Bin and Flight_category if present
extra_groups <- c("Time_Bin","Flight_category")
all_groups   <- unique(c(avail_groups, extra_groups[extra_groups %in% names(sp)]))

# ── C1: All distance matrices ─────────────────────────────────────────────────
# Build the distance tables used for the comparison tests.
# (a) shapePC1 + shapePC2 (already defined above)
shape_pc_cols <- grep("^shapePC", names(sp), value=TRUE)
cat(sprintf("  Shape PC axes available: %s\n", paste(shape_pc_cols, collapse=", ")))

shape_dist <- dist(sp[,shape_pc_cols], method="euclidean")

# (b) All shape PCs together (same as above if only PC1+PC2; includes more if available)
# Already covered by shape_dist

# (c) 6 biomechanical metrics (z-scored Euclidean)
bio_present <- BIO_VARS[BIO_VARS %in% names(sp)]
bio_scaled_mat <- scale(sp[,bio_present])
bio_dist    <- dist(bio_scaled_mat, method="euclidean")

# (d) WAR + r2_hat + von_mises_stress
war_vars    <- c("aspect_ratio","r2_hat","von_mises_stress")
war_vars    <- war_vars[war_vars %in% names(sp)]
war_dist    <- dist(scale(sp[,war_vars]), method="euclidean")

# (e) shape PCs + biomechanical metrics combined
shape_bio_mat  <- cbind(scale(sp[,shape_pc_cols]), bio_scaled_mat)
shape_bio_dist <- dist(shape_bio_mat, method="euclidean")

dist_list <- list(
  shapePC_only    = shape_dist,
  biomech_6vars   = bio_dist,
  WAR_2MA_SVM     = war_dist,
  shapePC_biomech = shape_bio_dist
)
