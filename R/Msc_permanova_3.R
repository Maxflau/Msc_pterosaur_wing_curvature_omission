# --------------------------------------------------------------------------------
# Check that all the objects this script needs already exist before continuing.
# --------------------------------------------------------------------------------
required_objects <- c("df", "pca_dist", "summarise_pcs", "test_dispersion",
                      "perm_clade", "perm_dep", "perm_env",
                      "perm_diet1", "perm_diet2", "perm_dietcombo", "perm_full",
                      "D_FULL", "MV_META")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop("Missing: ", paste(missing_objects, collapse = ", "),
       "\n  Source Msc_stat_permanova_2.R AND Msc_stat_permanova.R",
       "\n  (via Msc_stats_00_setup.R) before this script.")
}

# The groups compared in every test below.
GROUPS <- c("clade", "Depositional", "Palaeoenvironment",
            "Diet_primary", "Diet_secondary", "Diet_combined")

perm_results_pca <- list(
  clade              = perm_clade,
  Depositional       = perm_dep,
  Palaeoenvironment  = perm_env,
  Diet_primary       = perm_diet1,
  Diet_secondary     = perm_diet2,
  Diet_combined      = perm_dietcombo)

# Turn one PERMANOVA (adonis2) result into a tidy, labelled table row set.
tidy_adonis2 <- function(res, factor_label) {
  tbl           <- as.data.frame(res)
  tbl$Term      <- rownames(tbl)
  tbl$Factor    <- factor_label
  rownames(tbl) <- NULL
  tbl <- tbl[, c("Factor", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
  names(tbl)[names(tbl) == "Pr(>F)"] <- "p_value"
  tbl
}

# Turn one betadisper (variance homogeneity) test into a tidy, labelled table row set.
tidy_betadisper <- function(dist_mat, groups, factor_label) {
  bd            <- betadisper(dist_mat, groups)
  pv            <- permutest(bd, permutations = 9999)
  tbl           <- as.data.frame(pv$tab)
  tbl$Term      <- rownames(tbl)
  tbl$Factor    <- factor_label
  rownames(tbl) <- NULL
  tbl <- tbl[, c("Factor", "Term", "Df", "Sum Sq", "Mean Sq", "F", "N.Perm", "Pr(>F)")]
  names(tbl)[names(tbl) == "Pr(>F)"]  <- "p_value"
  names(tbl)[names(tbl) == "Sum Sq"]  <- "SumSq"
  names(tbl)[names(tbl) == "Mean Sq"] <- "MeanSq"
  tbl
}

OUT <- "output/results/"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("\n\n=== EXPORTING CSV FILES ===\n")

# --- 01  PC summaries by group ------------------------------------------------

for (g in GROUPS) {
  pc_tbl <- summarise_pcs(df, g)
  write.csv(pc_tbl,
            file      = paste0(OUT, "01_summary_pcs_", g, ".csv"),
            row.names = FALSE)
}
cat("01 — PC summary CSVs written.\n")


# --- 03  Individual PERMANOVA results — PCA space ----------------------------

pca_perm_list <- lapply(names(perm_results_pca), function(g) {
  tidy_adonis2(perm_results_pca[[g]], g)
})
pca_perm_df <- do.call(rbind, pca_perm_list)
write.csv(pca_perm_df,
          file      = paste0(OUT, "03_permanova_pca.csv"),
          row.names = FALSE)
cat("03 — PCA PERMANOVA CSV written.\n")


# --- 04  Full marginal model — PCA space -------------------------------------

write.csv(tidy_adonis2(perm_full, "Full model (PCA)"),
          file      = paste0(OUT, "04_permanova_pca_full_model.csv"),
          row.names = FALSE)
cat("04 — PCA full model CSV written.\n")


# --- 07  Betadisper — PCA space ----------------------------------------------

set.seed(42)

pca_disp_list <- lapply(GROUPS, function(g) tidy_betadisper(pca_dist, df[[g]], g))
pca_disp_df <- do.call(rbind, pca_disp_list)
write.csv(pca_disp_df,
          file      = paste0(OUT, "07_betadisper_pca.csv"),
          row.names = FALSE)
cat("07 — PCA betadisper CSV written.\n")


# --------------------------------------------------------------------------------
# 08. Betadisper - biomechanical (multivariate) space.
# --------------------------------------------------------------------------------
bio_dist <- D_FULL

bio_groups <- intersect(GROUPS, colnames(MV_META))
skipped_groups <- setdiff(GROUPS, bio_groups)
if (length(skipped_groups) > 0) {
  cat("NOTE: skipping betadisper (biomech) for groups not in MV_META:",
      paste(skipped_groups, collapse = ", "), "\n")
}

bio_disp_list <- lapply(bio_groups, function(g) {
  grp <- as.character(MV_META[[g]])
  ok <- !is.na(grp) & grp != ""
  if (sum(ok) < 10 || length(unique(grp[ok])) < 2) {
    cat("SKIPPED betadisper (biomech) for", g, "- too few usable specimens or levels\n")
    return(NULL)
  }
  Ds <- as.dist(as.matrix(bio_dist)[ok, ok])
  tidy_betadisper(Ds, factor(grp[ok]), g)
})
bio_disp_df <- do.call(rbind, bio_disp_list)

write.csv(bio_disp_df, file      = paste0(OUT, "08_betadisper_biomech.csv"), row.names = FALSE)
cat("08 — Biomechanical betadisper CSV written.\n")


# --------------------------------------------------------------------------------
# 09. Kruskal-Wallis tests, univariate, per metric x group.
# --------------------------------------------------------------------------------
KW_METRICS <- intersect(c("aspect_ratio", "r2_hat", "von_mises_stress",
                          "wing_loading_ratio", "wing_curvature", "shape_complexity",
                          "reynolds", "PC1", "PC2"), colnames(df))
cat("\nKruskal-Wallis metrics:", paste(KW_METRICS, collapse = ", "), "\n\n")

kw_results <- data.frame()
for (g in GROUPS) {
  for (v in KW_METRICS) {
    ok <- !is.na(df[[v]]) & !is.na(df[[g]]) & df[[g]] != ""
    if (sum(ok) < 10 || length(unique(df[[g]][ok])) < 2) next
    kw <- kruskal.test(df[[v]][ok], as.factor(df[[g]][ok]))
    kw_results <- rbind(kw_results, data.frame(
      Group = g, Variable = v,
      H = round(kw$statistic, 3), df_kw = kw$parameter,
      p_value = round(kw$p.value, 4),
      sig = ifelse(kw$p.value < 0.001, "***",
                   ifelse(kw$p.value < 0.01,  "**",
                          ifelse(kw$p.value < 0.05,  "*",
                                 ifelse(kw$p.value < 0.1,   ".", ""))))))
  }
}

cat("\nKruskal-Wallis omnibus results:\n")
print(kw_results, row.names = FALSE)

write.csv(kw_results, file      = paste0(OUT, "09_kruskal_wallis.csv"), row.names = FALSE)
cat("09 — Kruskal-Wallis CSV written.\n\n")
