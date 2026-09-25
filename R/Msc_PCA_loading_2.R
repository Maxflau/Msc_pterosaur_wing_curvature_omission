# Use the dplyr versions of common helper names.
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange

# --- 1. Load the PCA table and tidy the group names ---------------------------
df <- read.csv("output/results/Supplemental_performance/performance_metrics_complete_CLADE.csv", stringsAsFactors = FALSE) %>%
  rename(Depositional = Depositional.settings.paleoenvironment,
         Palaeoenvironment = Environment, Diet_primary = Diet.1, Diet_secondary = Diet.2)

# Build the combined diet label if it is missing from the table.
if (!"Diet_combined" %in% names(df) && all(c("Diet_primary", "Diet_secondary") %in% names(df))) {
  d1 <- trimws(as.character(df$Diet_primary)); d2 <- trimws(as.character(df$Diet_secondary))
  df$Diet_combined <- ifelse(is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = " + "))
}

# second_moment renamed to r2_hat (non-dimensional second moment of area) -
# wing_loading_ratio added (size-corrected, positive; replaces wing_loading
# in the multivariate PCA itself since Msc_impossible_regions_v3.R).
# Choose the biomechanical columns and group columns used below.
BIO_VARS <- c("aspect_ratio","r2_hat","wing_loading_ratio",
              "von_mises_stress","wing_curvature","shape_complexity")
GROUPS   <- c("clade","Depositional","Palaeoenvironment",
              "Diet_primary","Diet_secondary","Diet_combined")

# --- 2. VARIANCE EXPLAINED ----------------------------------------------------
# Read from the ACTUAL current PCA (pca_performance, from
# Msc_impossible_regions_v3.R) rather than hard-coded percentages - the old
# 47.70%/20.94% values were computed on a different metric set (before
# wing_loading_ratio replaced wing_loading, before tip_angle/pitch_agility
# were added then removed again) and no longer describe this PCA.
if (exists("pca_performance")) {
  imp <- summary(pca_performance)$importance
  var_exp <- c(PC1 = round(imp[2, 1] * 100, 2), PC2 = round(imp[2, 2] * 100, 2))
  cum_var <- c(PC1 = var_exp[["PC1"]], PC2 = round(sum(imp[2, 1:2]) * 100, 2))
} else {
  stop("pca_performance not found - source Msc_impossible_regions_v3.R first ",
       "so the variance explained is read from the CURRENT PCA, not a stale hard-coded value.")
}
var_tbl <- data.frame(PC = names(var_exp), Var_pct = var_exp, Cumul_pct = cum_var)
cat("=== VARIANCE EXPLAINED (PC1 + PC2) ===\n"); print(var_tbl)
write.csv(var_tbl, "output/results/Supplemental_performance/01_pca_variance_explained.csv", row.names = FALSE)

# --- 2B. Recalculate how strongly each variable lines up with the PCs ---------
# Scale the selected measurements so they are on a comparable scale.
bio_scaled_mat <- scale(df[, BIO_VARS])

loadings_computed <- data.frame(
  variable = BIO_VARS,
  PC1      = apply(bio_scaled_mat, 2, function(x) cor(x, df$PC1, use = "complete.obs")),
  PC2      = apply(bio_scaled_mat, 2, function(x) cor(x, df$PC2, use = "complete.obs"))
)

cat("\n=== RECALCULATED LOADINGS (correlation with PC scores) ===\n")
print(loadings_computed)
write.csv(loadings_computed, "output/results/Supplemental_performance/00_pca_loadings_recalculated.csv", row.names = FALSE)

load_df <- loadings_computed

# --- 3. LOADING BIPLOT (PC1 vs PC2) ------------------------------------------
arrow_scale <- max(abs(df[, c("PC1","PC2")])) * 0.6

p_biplot <- ggplot() +
  geom_point(data = df, aes(PC1, PC2, colour = clade), size = 1.5, alpha = 0.7) +
  geom_segment(data = load_df,
               aes(x=0, y=0, xend=PC1*arrow_scale, yend=PC2*arrow_scale),
               arrow = arrow(length = unit(0.2,"cm")), colour="black", linewidth=0.6) +
  geom_label_repel(data = load_df,
                   aes(x=PC1*arrow_scale, y=PC2*arrow_scale, label=variable),
                   size=3, max.overlaps=20) +
  labs(title  = "Biomechanical morphospace — PC1 vs PC2",
       x = paste0("PC1 (",var_exp[["PC1"]],"%)"), y = paste0("PC2 (",var_exp[["PC2"]],"%)"),
       colour = "Clade") + theme_bw()
ggsave("output/plots_PDF/supplementals/02_biplot_PC1_PC2.pdf", p_biplot, width=12, height=8)

# --- 4. Summarise PC1 and PC2 within each group -------------------------------
summarise_pcs <- function(data, g) {
  data %>% group_by(across(all_of(g))) %>%
    summarise(N=n(), PC1_mean=round(mean(PC1),3), PC1_sd=round(sd(PC1),3),
              PC2_mean=round(mean(PC2),3), PC2_sd=round(sd(PC2),3), .groups="drop")
}
for (g in GROUPS) {
  tbl <- summarise_pcs(df, g)
  cat("\n--- PC summary:", g, "---\n")
  # Show a fixed number of decimals so every row lines up.
  tbl_display <- tbl
  tbl_display$PC1_mean <- format_fixed(tbl$PC1_mean, 3)
  tbl_display$PC1_sd   <- format_fixed(tbl$PC1_sd, 3)
  tbl_display$PC2_mean <- format_fixed(tbl$PC2_mean, 3)
  tbl_display$PC2_sd   <- format_fixed(tbl$PC2_sd, 3)
  print(tbl_display)
  write.csv(tbl_display, paste0("output/results/Supplemental_performance/03_summary_pcs_", g, ".csv"), row.names = FALSE)
}

# --- 5. PERMANOVA ON PC1 + PC2 -----------------------------------------------
pc_dist <- dist(df[, c("PC1","PC2")], method = "euclidean")
set.seed(42)

perm_rows <- lapply(GROUPS, function(g) {
  res <- adonis2(as.formula(paste("pc_dist ~", g)),
                 data=df, permutations=9999, method="euclidean")
  tbl <- as.data.frame(res); tbl$Term <- rownames(tbl); tbl$Factor <- g
  tbl[!is.na(tbl$F), ]
})
perm_df <- do.call(rbind, perm_rows) %>%
  select(Factor, Term, Df, SumOfSqs, R2, F, `Pr(>F)`) %>%
  rename(p_value=`Pr(>F)`) %>%
  mutate(sig = case_when(p_value<0.001~"***", p_value<0.01~"**",
                         p_value<0.05~"*",   p_value<0.1~".", TRUE~""))
cat("\n=== PERMANOVA (PC1+PC2) ===\n")
# Show a fixed number of decimals so every row lines up.
perm_df_display <- perm_df
perm_df_display$SumOfSqs <- format_fixed(perm_df$SumOfSqs, 4)
perm_df_display$R2       <- format_fixed(perm_df$R2, 4)
perm_df_display$F        <- format_fixed(perm_df$F, 3)
perm_df_display$p_value  <- format_fixed(perm_df$p_value, 4)
print(perm_df_display)
write.csv(perm_df_display, "output/results/Supplemental_performance/04_permanova_pc1pc2.csv", row.names = FALSE)

# --- 6. PAIRWISE PERMANOVA — BONFERRONI (clades & diet combo) ----------------
set.seed(42)
for (g in c("clade","Diet_combined")) {
  pw <- pairwise.adonis2(as.formula(paste("pc_dist ~", g)),
                         data=df, permutations=9999, p.adjust.m="bonferroni")
  pw_rows <- lapply(names(pw)[names(pw) != "parent_call"], function(nm) {
    tbl <- as.data.frame(pw[[nm]]); tbl$Comparison <- nm; tbl$Term <- rownames(tbl); tbl
  })
  pw_df <- do.call(rbind, pw_rows) %>%
    filter(!Term %in% c("Residual","Total")) %>%
    select(Comparison, Df, SumOfSqs, R2, F, `Pr(>F)`) %>%
    rename(p_adj=`Pr(>F)`) %>%
    mutate(sig = case_when(p_adj<0.001~"***", p_adj<0.01~"**",
                           p_adj<0.05~"*",   p_adj<0.1~".", TRUE~""))
  cat("\n--- Pairwise PERMANOVA:", g, "---\n")
  # Show a fixed number of decimals so every row lines up.
  pw_df_display <- as.data.frame(pw_df)
  pw_df_display$SumOfSqs <- format_fixed(pw_df$SumOfSqs, 4)
  pw_df_display$R2       <- format_fixed(pw_df$R2, 4)
  pw_df_display$F        <- format_fixed(pw_df$F, 3)
  pw_df_display$p_adj    <- format_fixed(pw_df$p_adj, 4)
  print(pw_df_display)
  write.csv(pw_df_display, paste0("output/results/Supplemental_performance/05_pairwise_permanova_", g, ".csv"), row.names=FALSE)
}

# --- 7. HOMOGENEITY OF DISPERSION --------------------------------------------
disp_rows <- lapply(GROUPS, function(g) {
  bd  <- betadisper(pc_dist, df[[g]])
  pv  <- permutest(bd, permutations=9999)
  tbl <- as.data.frame(pv$tab); tbl$Term <- rownames(tbl); tbl$Factor <- g
  tbl[!is.na(tbl$F), ]
})
disp_df <- do.call(rbind, disp_rows) %>%
  rename(p_value=`Pr(>F)`) %>%
  select(Factor, Term, Df, `Sum Sq`, `Mean Sq`, F, p_value)
# Show a fixed number of decimals so every row lines up.
disp_df_display <- disp_df
disp_df_display$`Sum Sq`  <- format_fixed(disp_df$`Sum Sq`, 4)
disp_df_display$`Mean Sq` <- format_fixed(disp_df$`Mean Sq`, 4)
disp_df_display$F         <- format_fixed(disp_df$F, 3)
disp_df_display$p_value   <- format_fixed(disp_df$p_value, 4)
write.csv(disp_df_display, "output/results/Supplemental_performance/06_betadisper_pc1pc2.csv", row.names=FALSE)

cat("\n=== All CSV files saved. Analysis complete. ===\n")
