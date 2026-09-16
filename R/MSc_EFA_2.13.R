# ── C2: PERMANOVA — all groups × all matrices ─────────────────────────────────
set.seed(42)
all_perm_rows <- list()

for (dname in names(dist_list)) {
  d <- dist_list[[dname]]
  cat(sprintf("\n  PERMANOVA on %s:\n", dname))
  for (g in all_groups) {
    if (!g %in% names(sp)) next
    g_clean <- sp[[g]]; g_clean[is.na(g_clean) | g_clean==""] <- NA
    if (sum(!is.na(g_clean)) < 10 || length(unique(na.omit(g_clean))) < 2) next
    sp_tmp <- sp; sp_tmp[[g]] <- g_clean
    sp_tmp  <- sp_tmp[!is.na(g_clean),]
    d_sub   <- as.dist(as.matrix(d)[!is.na(g_clean), !is.na(g_clean)])
    res <- tryCatch(
      adonis2(as.formula(paste("d_sub ~", g)), data=sp_tmp,
              permutations=9999, method="euclidean"),
      error=function(e) NULL
    )
    if (!is.null(res)) {
      row <- tidy_perm(res, g); row$Matrix <- dname
      all_perm_rows[[paste(dname, g)]] <- row
      cat(sprintf("    %s: R2=%.3f F=%.1f p=%s\n",
                  g, row$R2[1], row$F[1],
                  ifelse(row$p_value[1]<0.001,"<0.001",round(row$p_value[1],3))))
    }
  }
}
all_perm_df <- do.call(rbind, all_perm_rows)
save_csv(all_perm_df, "EFA_permanova_all_matrices_all_groups")

# Subset tables per matrix for readability
for (dname in names(dist_list)) {
  sub <- all_perm_df[all_perm_df$Matrix==dname, ]
  save_csv(sub, paste0("EFA_permanova_", dname))
}

# ── C3: Full marginal model (shape PCs only) ─────────────────────────────────
marg_vars <- c("clade","Depositional","Palaeoenvironment","Diet_combo")
marg_vars <- marg_vars[marg_vars %in% names(sp)]
if (length(marg_vars) >= 2) {
  perm_full <- adonis2(
    as.formula(paste("shape_dist ~", paste(marg_vars, collapse="+"))),
    data=sp, permutations=9999, method="euclidean", by="margin"
  )
  save_csv(tidy_perm(perm_full,"Full marginal (shapePC)"),
           "EFA_permanova_full_marginal_shapePC")
}

# ── C4: PAIRWISE PERMANOVA — all requested groups × shape_dist ───────────────
# Groups: clade, Diet_primary, Diet_secondary, Palaeoenvironment,
#         Depositional, Time_Bin, Diet_combo
pairwise_groups <- c("clade","Diet_primary","Diet_secondary",
                     "Palaeoenvironment","Depositional","Time_Bin","Diet_combo")
pairwise_groups <- pairwise_groups[pairwise_groups %in% names(sp)]

cat("\n  Pairwise PERMANOVA (Bonferroni) on shapePC space:\n")
all_pw_list <- list()

for (g in pairwise_groups) {
  g_vals <- sp[[g]]; g_vals[is.na(g_vals)|g_vals==""] <- NA
  n_levels <- length(unique(na.omit(g_vals)))
  if (n_levels < 2) { cat(sprintf("    Skipping %s: <2 levels\n", g)); next }
  if (n_levels > 30) {
    cat(sprintf("    Warning: %s has %d levels — pairwise may be slow\n", g, n_levels))
  }
  sp_pw <- sp[!is.na(g_vals),]
  d_pw  <- as.dist(as.matrix(shape_dist)[!is.na(g_vals), !is.na(g_vals)])
  
  pw <- tryCatch(
    pairwise.adonis2(as.formula(paste("d_pw ~", g)), data=sp_pw,
                     permutations=9999, p.adjust.m="bonferroni"),
    error=function(e) { cat(sprintf("    %s failed: %s\n",g,e$message)); NULL }
  )
  if (!is.null(pw)) {
    pw_df <- tidy_pairwise(pw, g)
    all_pw_list[[g]] <- pw_df
    save_csv(as.data.frame(pw_df), paste0("EFA_pairwise_permanova_", g))
    cat(sprintf("    %s: %d comparisons\n", g, nrow(pw_df)))
  }
}

# Summary: significant pairwise comparisons across all groups
if (length(all_pw_list) > 0) {
  all_pw_df <- do.call(rbind, all_pw_list)
  save_csv(as.data.frame(all_pw_df), "EFA_pairwise_permanova_ALL_groups")
  cat(sprintf("\n  Total pairwise comparisons: %d | Significant (p<0.05): %d\n",
              nrow(all_pw_df), sum(all_pw_df$p_adj<0.05, na.rm=TRUE)))
}


# =============================================================================
# SECTION D — BETADISPER: HOMOGENEITY OF DISPERSION (Msc_Perma_3)
# =============================================================================
cat("\n── Section D: Betadisper ──\n")

disp_rows <- list()
for (dname in names(dist_list)) {
  for (g in all_groups) {
    if (!g %in% names(sp)) next
    g_clean <- sp[[g]]; g_clean[is.na(g_clean)|g_clean==""] <- NA
    if (sum(!is.na(g_clean)) < 6 || length(unique(na.omit(g_clean))) < 2) next
    d_sub <- as.dist(as.matrix(dist_list[[dname]])[!is.na(g_clean), !is.na(g_clean)])
    bd <- tryCatch(betadisper(d_sub, na.omit(g_clean)), error=function(e) NULL)
    if (is.null(bd)) next
    pv <- tryCatch(permutest(bd, permutations=999), error=function(e) NULL)
    if (is.null(pv)) next
    tbl <- as.data.frame(pv$tab); tbl$Term <- rownames(tbl)
    tbl$Factor <- g; tbl$Matrix <- dname
    disp_rows[[paste(dname,g)]] <- tbl[!is.na(tbl$F),]
  }
}
disp_df <- do.call(rbind, disp_rows)
save_csv(disp_df, "EFA_betadisper_all_matrices")