# ── 9A. DIRECTORIES ───────────────────────────────────────────────────────────
# Make sure the output folders exist before saving results.
res_dir  <- "output/results/Supplemental_EFA"
stat_dir <- res_dir                       # save_csv() writes here
plot_dir <- "output/plots_PDF"
png_dir  <- "output/plots"
for (dd in c(res_dir, plot_dir, png_dir)) dir.create(dd, showWarnings=FALSE, recursive=TRUE)
cat("Output directories ready\n")

# ── 9B. DATA PREPARATION ──────────────────────────────────────────────────────
# Prepare one clean table with the names and groups used later.
stat_df <- shape_perf

rename_safe <- function(df, old, new) {
  if (old %in% names(df) && !new %in% names(df)) dplyr::rename(df, !!new := !!sym(old)) else df
}
stat_df <- stat_df %>%
  rename_safe("Depositional.settings.paleoenvironment", "Depositional") %>%
  rename_safe("Environment",                            "Palaeoenvironment") %>%
  rename_safe("Diet.1",                                 "Diet_primary") %>%
  rename_safe("Diet.2",                                 "Diet_secondary")

# Diet_combo is built HERE, on stat_df, after the renames. Building it on
# shape_perf earlier does not help: stat_df is a copy taken before that point,
# which is what produced "object 'Diet_combo' not found" in adonis2.
if (!"Diet_combo" %in% names(stat_df) && "Diet_primary" %in% names(stat_df)) {
  d1 <- trimws(as.character(stat_df$Diet_primary))
  d2 <- if ("Diet_secondary" %in% names(stat_df))
    trimws(as.character(stat_df$Diet_secondary)) else rep(NA_character_, length(d1))

  # The pair is sorted so "A+B" and "B+A" form one category. The raw columns
  # contain both orderings, which otherwise splits the counts across duplicates.
  stat_df$Diet_combo <- mapply(function(a, b) {
    if (is.na(a) || a == "") return(NA_character_)
    if (is.na(b) || b == "" || b == a) a else paste(sort(c(a, b)), collapse = "+")
  }, d1, d2, USE.NAMES = FALSE)

  tb <- table(stat_df$Diet_combo)
  cat(sprintf("Diet_combo: %d levels, %d with n >= 3\n", length(tb), sum(tb >= 3)))
}

# Time bins, needed by Section 10
if (!"Time_Bin" %in% names(stat_df)) {
  pcol <- intersect(c("Period.Name", "Period_Name"), names(stat_df))[1]
  if (!is.na(pcol)) {
    pc <- trimws(as.character(stat_df[[pcol]]))
    pc[pc %in% c("Early Jurassic", "Middle Jurassic")] <- "Early+Middle Jurassic"
    stat_df$Time_Bin <- factor(pc,
                               levels = c("Late Triassic","Early+Middle Jurassic","Late Jurassic",
                                          "Early Cretaceous","Late Cretaceous"))
    shape_perf$Time_Bin <- stat_df$Time_Bin   # Section 10 reads shape_perf
    cat(sprintf("Time_Bin: %d specimens, %d bins\n",
                sum(!is.na(stat_df$Time_Bin)), nlevels(droplevels(stat_df$Time_Bin))))
  }
}

BIO_VARS <- c("aspect_ratio","wing_loading","von_mises_stress",
              "wing_curvature","shape_complexity","r2_hat","pareto_rank_ratio")
BIO_VARS <- BIO_VARS[BIO_VARS %in% names(stat_df)]

# Only factors actually present, so no formula can reference a missing column
GROUPS <- c("clade","Depositional","Palaeoenvironment",
            "Diet_primary","Diet_secondary","Diet_combo")
GROUPS <- GROUPS[GROUPS %in% names(stat_df)]

cat("Metrics:", paste(BIO_VARS, collapse=", "), "\n")
cat("Factors:", paste(GROUPS, collapse=", "), "\n")

# ── 9C. SUMMARY STATISTICS ────────────────────────────────────────────────────
# Save simple summaries for the main groups.
summarise_group <- function(data, grp_col) {
  data %>%
    filter(!is.na(.data[[grp_col]]), .data[[grp_col]] != "") %>%
    group_by(.data[[grp_col]]) %>%
    summarise(n = n(),
              across(all_of(BIO_VARS),
                     list(mean = ~round(mean(.x, na.rm=TRUE), 4),
                          sd   = ~round(sd(.x,   na.rm=TRUE), 4)),
                     .names = "{.col}__{.fn}"),
              shapePC1_mean = round(mean(shapePC1, na.rm=TRUE), 4),
              shapePC1_sd   = round(sd(shapePC1,   na.rm=TRUE), 4),
              shapePC2_mean = round(mean(shapePC2, na.rm=TRUE), 4),
              shapePC2_sd   = round(sd(shapePC2,   na.rm=TRUE), 4),
              .groups = "drop") %>%
    rename(group = 1)
}

for (g in intersect(c("clade","Order","Diet_combo"), names(stat_df))) {
  write.csv(summarise_group(stat_df, g),
            file.path(res_dir, paste0("summary_statistics_by_", tolower(g), ".csv")),
            row.names = FALSE)
  cat("Written: summary_statistics_by_", tolower(g), ".csv\n", sep="")
}

# ── 9D. PERMANOVA ON SHAPE PC SPACE ───────────────────────────────────────────
# Test whether groups differ in EFA shape space.
shape_dist <- dist(as.matrix(stat_df[, c("shapePC1","shapePC2")]), method="euclidean")
set.seed(42)

tidy_perm <- function(res, factor_label) {
  tbl <- as.data.frame(res); tbl$Term <- rownames(tbl); tbl$Factor <- factor_label
  tbl <- tbl[!is.na(tbl$F), c("Factor","Term","Df","SumOfSqs","R2","F","Pr(>F)")]
  names(tbl)[names(tbl)=="Pr(>F)"] <- "p_value"
  tbl
}

perm_shape_df <- do.call(rbind, Filter(Negate(is.null), lapply(GROUPS, function(g) {
  sub <- stat_df[!is.na(stat_df[[g]]) & stat_df[[g]] != "", ]
  if (length(unique(sub[[g]])) < 2) return(NULL)
  d <- dist(as.matrix(sub[, c("shapePC1","shapePC2")]), method="euclidean")
  tidy_perm(adonis2(as.formula(paste("d ~", g)), data=sub, permutations=9999), g)
}))) %>%
  mutate(sig = case_when(p_value<0.001~"***", p_value<0.01~"**",
                         p_value<0.05~"*", p_value<0.1~".", TRUE~""))

write.csv(perm_shape_df, file.path(res_dir,"permanova_shape_pcs.csv"), row.names=FALSE)
cat("Written: permanova_shape_pcs.csv\n")

# Marginal model, built only from factors present and complete. adonis2 drops no
# rows itself, so any NA in a term would silently unbalance the design.
MARGIN_TERMS <- intersect(c("clade","Depositional","Palaeoenvironment","Diet_combo"),
                          GROUPS)
if (length(MARGIN_TERMS) >= 2) {
  cc <- complete.cases(stat_df[, MARGIN_TERMS])
  sub <- stat_df[cc, ]
  cat(sprintf("Marginal model on %d of %d specimens: %s\n",
              sum(cc), nrow(stat_df), paste(MARGIN_TERMS, collapse=" + ")))
  d <- dist(as.matrix(sub[, c("shapePC1","shapePC2")]), method="euclidean")
  perm_shape_full <- adonis2(as.formula(paste("d ~", paste(MARGIN_TERMS, collapse=" + "))),
                             data=sub, permutations=9999, by="margin")
  write.csv(tidy_perm(perm_shape_full, "marginal"),
            file.path(res_dir,"permanova_shape_marginal.csv"), row.names=FALSE)
  cat("Written: permanova_shape_marginal.csv\n")
} else {
  cat("Fewer than two factors available - marginal model skipped.\n")
}
