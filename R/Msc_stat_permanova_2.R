df <- read.csv("clade_analyses_v4/results/performance_metrics_complete_CLADE.csv",
               stringsAsFactors = FALSE)

# Rename long column for convenience
df <- df %>%
  rename(Depositional     = Depositional.settings.paleoenvironment,
         Palaeoenvironment = Environment,
         Diet_primary     = Diet.1,
         Diet_secondary   = Diet.2)

# Diet_combined is built in Msc_stat_prep.R and should already be in the CSV
# if that script ran before the export - built here defensively otherwise.
if (!"Diet_combined" %in% names(df)) {
  d1 <- trimws(as.character(df$Diet_primary)); d2 <- trimws(as.character(df$Diet_secondary))
  df$Diet_combined <- ifelse(is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = " + "))
}

# Quick check
cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("NAs in PC1:", sum(is.na(df$PC1)),
    "| NAs in PC2:", sum(is.na(df$PC2)), "\n\n")


# --- 2. PC SCORE SUMMARY BY GROUP --------------------------------------------
summarise_pcs <- function(data, group_var) {
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      N       = n(),
      PC1_mean   = round(mean(PC1),  3),
      PC1_sd     = round(sd(PC1),    3),
      PC1_median = round(median(PC1),3),
      PC2_mean   = round(mean(PC2),  3),
      PC2_sd     = round(sd(PC2),    3),
      PC2_median = round(median(PC2),3),
      .groups = "drop"
    )
}

cat("=== PC1 & PC2 SUMMARY BY CLADE ===\n")
sum_clade <- summarise_pcs(df, "clade")
print(kable(sum_clade, format = "simple"))

cat("\n=== PC1 & PC2 SUMMARY BY DEPOSITIONAL SETTING ===\n")
sum_dep <- summarise_pcs(df, "Depositional")
print(kable(sum_dep, format = "simple"))

cat("\n=== PC1 & PC2 SUMMARY BY PALAEOENVIRONMENT ===\n")
sum_env <- summarise_pcs(df, "Palaeoenvironment")
print(kable(sum_env, format = "simple"))

cat("\n=== PC1 & PC2 SUMMARY BY PRIMARY DIET (Diet.1) ===\n")
sum_diet1 <- summarise_pcs(df, "Diet_primary")
print(kable(sum_diet1, format = "simple"))

cat("\n=== PC1 & PC2 SUMMARY BY SECONDARY DIET (Diet.2) ===\n")
sum_diet2 <- summarise_pcs(df, "Diet_secondary")
print(kable(sum_diet2, format = "simple"))

cat("\n=== PC1 & PC2 SUMMARY BY DIET COMBINATION ===\n")
sum_dietcombo <- summarise_pcs(df, "Diet_combined")
print(kable(sum_dietcombo, format = "simple"))


# --- 3. BUILD DISTANCE MATRIX FOR PERMANOVA ----------------------------------
pca_mat  <- as.matrix(df[, c("PC1", "PC2")])
pca_dist <- dist(pca_mat, method = "euclidean")


# --- 4. PERMANOVA (adonis2) --------------------------------------------------
set.seed(42)

cat("\n\n=== PERMANOVA: PCA SPACE ~ CLADE ===\n")
perm_clade <- adonis2(pca_dist ~ clade,
                      data         = df,
                      permutations = 9999,
                      method       = "euclidean")
print(perm_clade)

cat("\n=== PERMANOVA: PCA SPACE ~ DEPOSITIONAL SETTING ===\n")
perm_dep <- adonis2(pca_dist ~ Depositional,
                    data         = df,
                    permutations = 9999,
                    method       = "euclidean")
print(perm_dep)

cat("\n=== PERMANOVA: PCA SPACE ~ PALAEOENVIRONMENT ===\n")
perm_env <- adonis2(pca_dist ~ Palaeoenvironment,
                    data         = df,
                    permutations = 9999,
                    method       = "euclidean")
print(perm_env)

cat("\n=== PERMANOVA: PCA SPACE ~ PRIMARY DIET (Diet.1) ===\n")
perm_diet1 <- adonis2(pca_dist ~ Diet_primary,
                      data         = df,
                      permutations = 9999,
                      method       = "euclidean")
print(perm_diet1)

cat("\n=== PERMANOVA: PCA SPACE ~ SECONDARY DIET (Diet.2) ===\n")
perm_diet2 <- adonis2(pca_dist ~ Diet_secondary,
                      data         = df,
                      permutations = 9999,
                      method       = "euclidean")
print(perm_diet2)

cat("\n=== PERMANOVA: PCA SPACE ~ DIET COMBINATION ===\n")
perm_dietcombo <- adonis2(pca_dist ~ Diet_combined,
                          data         = df,
                          permutations = 9999,
                          method       = "euclidean")
print(perm_dietcombo)

cat("\n=== PERMANOVA: FULL MODEL (Clade + Depositional + Palaeoenvironment + Diet combination) ===\n")
perm_full <- adonis2(pca_dist ~ clade + Depositional + Palaeoenvironment + Diet_combined,
                     data         = df,
                     permutations = 9999,
                     method       = "euclidean",
                     by           = "margin")
print(perm_full)


# --- 5. HOMOGENEITY OF DISPERSION (betadisper) -------------------------------
test_dispersion <- function(dist_mat, groups, label) {
  bd  <- betadisper(dist_mat, groups)
  pv  <- permutest(bd, permutations = 9999)
  cat("\n--- Dispersion test:", label, "---\n")
  print(pv)
}

cat("\n\n=== HOMOGENEITY OF DISPERSION TESTS ===\n")
test_dispersion(pca_dist, df$clade,             "Clade")
test_dispersion(pca_dist, df$Depositional,      "Depositional setting")
test_dispersion(pca_dist, df$Palaeoenvironment, "Palaeoenvironment")
test_dispersion(pca_dist, df$Diet_primary,      "Primary diet (Diet.1)")
test_dispersion(pca_dist, df$Diet_secondary,    "Secondary diet (Diet.2)")
test_dispersion(pca_dist, df$Diet_combined,     "Diet combination")


# --- 6. OPTIONAL: PAIRWISE PERMANOVA FOR CLADES ------------------------------
library(pairwiseAdonis)
cat("\n=== PAIRWISE PERMANOVA: CLADES ===\n")
pw_clade <- pairwise.adonis2(pca_dist ~ clade,
                             data         = df,
                             permutations = 9999)
print(pw_clade)