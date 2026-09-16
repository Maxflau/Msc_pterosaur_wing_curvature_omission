# Turn the PERMANOVA output into a simple table that can be saved.
tidy_adonis2 <- function(res, factor_label) {
  tbl           <- as.data.frame(res)
  tbl$Term      <- rownames(tbl)
  tbl$Factor    <- factor_label
  rownames(tbl) <- NULL
  tbl <- tbl[, c("Factor", "Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
  names(tbl)[names(tbl) == "Pr(>F)"] <- "p_value"
  tbl
}

# Turn the dispersion test output into a simple table that can be saved.
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

cat("\n\n=== EXPORTING CSV FILES ===\n")

# --- 01 & 02  PC and biomechanical summaries ---------------------------------
# Save the basic summary tables for each grouping variable.

for (g in GROUPS) {
  pc_tbl <- summarise_pcs(df, g)
  write.csv(pc_tbl,
            file      = paste0(OUT, "01_summary_pcs_", g, ".csv"),
            row.names = FALSE)

  bio_tbl <- bio_summary(df, g)
  write.csv(bio_tbl,
            file      = paste0(OUT, "02_summary_biomech_", g, ".csv"),
            row.names = FALSE)
}
cat("01 & 02 — Summary CSVs written.\n")

# --- 03  Individual PERMANOVA results — PCA space ----------------------------
# Test each grouping variable one at a time in PCA space.

set.seed(42)

pca_perm_list <- list()
for (g in GROUPS) {
  res <- adonis2(as.formula(paste("pca_dist ~", g)),
                 data         = df,
                 permutations = 9999,
                 method       = "euclidean")
  pca_perm_list[[g]] <- tidy_adonis2(res, g)
}

pca_perm_df <- do.call(rbind, pca_perm_list)
write.csv(pca_perm_df,
          file      = paste0(OUT, "03_permanova_pca.csv"),
          row.names = FALSE)
cat("03 — PCA PERMANOVA CSV written.\n")

# --- 04  Full marginal model — PCA space -------------------------------------
# Test all main grouping variables together in PCA space.

perm_pca_full_export <- adonis2(
  pca_dist ~ clade + Depositional + Palaeoenvironment + Diet_combined,
  data         = df,
  permutations = 9999,
  method       = "euclidean",
  by           = "margin"
)
write.csv(tidy_adonis2(perm_pca_full_export, "Full model (PCA)"),
          file      = paste0(OUT, "04_permanova_pca_full_model.csv"),
          row.names = FALSE)
cat("04 — PCA full model CSV written.\n")

# --- 05  Individual PERMANOVA results — biomechanical space ------------------
# Repeat the one-factor tests in biomechanical space.

bio_perm_list <- list()
for (g in GROUPS) {
  res <- adonis2(as.formula(paste("bio_dist ~", g)),
                 data         = df,
                 permutations = 9999,
                 method       = "euclidean")
  bio_perm_list[[g]] <- tidy_adonis2(res, g)
}

bio_perm_df <- do.call(rbind, bio_perm_list)
write.csv(bio_perm_df,
          file      = paste0(OUT, "05_permanova_biomech.csv"),
          row.names = FALSE)
cat("05 — Biomechanical PERMANOVA CSV written.\n")

# --- 06  Full marginal model — biomechanical space ---------------------------
# Test all main grouping variables together in biomechanical space.

perm_bio_full_export <- adonis2(
  bio_dist ~ clade + Depositional + Palaeoenvironment + Diet_combined,
  data         = df,
  permutations = 9999,
  method       = "euclidean",
  by           = "margin"
)
write.csv(tidy_adonis2(perm_bio_full_export, "Full model (Biomech)"), file  = paste0(OUT, "06_permanova_biomech_full_model.csv"),row.names = FALSE)
cat("06 — Biomechanical full model CSV written.\n")

# --- 07  Betadisper — PCA space ----------------------------------------------
# Check whether group spread differs in PCA space.

set.seed(42)
pca_disp_list <- list()
for (g in GROUPS) {
  pca_disp_list[[g]] <- tidy_betadisper(pca_dist, df[[g]], g)
}
pca_disp_df <- do.call(rbind, pca_disp_list)
write.csv(pca_disp_df, file= paste0(OUT, "07_betadisper_pca.csv"),row.names = FALSE)
cat("07 — PCA betadisper CSV written.\n")

# --- 08  Betadisper — biomechanical space ------------------------------------
# Check whether group spread differs in biomechanical space.
bio_disp_list <- list()
for (g in GROUPS) {
  bio_disp_list[[g]] <- tidy_betadisper(bio_dist, df[[g]], g)
}

bio_disp_df <- do.call(rbind, bio_disp_list)
write.csv(bio_disp_df, file= paste0(OUT, "08_betadisper_biomech.csv"),row.names = FALSE)
cat("08 — Biomechanical betadisper CSV written.\n")

# --- 09  Kruskal-Wallis results (already computed in section 10) -------------
# Save the non-parametric test results that were already created earlier.

write.csv(kw_results, file = paste0(OUT, "09_kruskal_wallis.csv"),row.names = FALSE)
cat("09 — Kruskal-Wallis CSV written.\n")
