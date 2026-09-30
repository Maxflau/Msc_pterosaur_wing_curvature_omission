# ── HELPERS ───────────────────────────────────────────────────────────────────
# Save files and set shared names used in the next tables.
save_csv  <- function(df, name) {
  path <- file.path(stat_dir, paste0(name, ".csv"))
  write.csv(df, path, row.names=FALSE)
  cat(sprintf("  ✓ %s.csv\n", name))
}
save_plot <- function(p, name, w=10, h=7) {
  pdf_path <- file.path(plot_dir, paste0(name, ".pdf"))
  png_path <- file.path(plot_dir, paste0(name, ".png"))
  ggsave(pdf_path, p, width=w, height=h, dpi=300)
  ggsave(png_path, p, width=w, height=h, dpi=150)
  cat(sprintf("  ✓ %s (PDF + PNG)\n", name))
}

# EFA axis labels
efa_lab1 <- if (exists("lab1")) lab1 else "EFA PC1"
efa_lab2 <- if (exists("lab2")) lab2 else "EFA PC2"

BIO_VARS <- c("aspect_ratio","r2_hat","wing_loading",
              "von_mises_stress","wing_curvature","shape_complexity")
GROUPS   <- c("clade","Depositional","Palaeoenvironment",
              "Diet_primary","Diet_secondary","Diet_combo")

# Use shape_perf if available; otherwise try loading from CSV
if (!exists("shape_perf")) {
  cat("shape_perf not found — attempting to load from CSV\n")
  sp_path <- file.path("EFA_stat-results",
                       "../output/results/specimen_biomechanical_pca_table.csv")
  if (file.exists(sp_path)) {
    shape_perf <- read.csv(sp_path, stringsAsFactors=FALSE)
  } else {
    stop("shape_perf not found. Source Msc_perf_surfaces_v4.R first.")
  }
}

# Rename environmental columns if present
ren_safe <- function(df, old, new) {
  if (old %in% names(df) && !new %in% names(df)) names(df)[names(df)==old] <- new; df
}
shape_perf <- shape_perf %>%
  ren_safe("Environment",                             "Palaeoenvironment") %>%
  ren_safe("Depositional.settings.paleoenvironment",  "Depositional") %>%
  ren_safe("Diet.1",                                  "Diet_primary") %>%
  ren_safe("Diet.2",                                  "Diet_secondary") %>%
  ren_safe("Period.Name",                             "Period_Name")

sp   <- shape_perf
cat(sprintf("Dataset: %d specimens | %d columns\n", nrow(sp), ncol(sp)))

# =============================================================================
# SECTION A — SUMMARY STATISTICS BY CLADE AND ORDER (Msc_Clade_sum.R)
# =============================================================================
# Make summary tables for clades, orders, and other key groups.
cat("\n── Section A: Summary statistics ──\n")

summarise_group <- function(data, grp_col) {
  data %>%
    filter(!is.na(.data[[grp_col]])) %>%
    group_by(.data[[grp_col]]) %>%
    summarise(
      n            = n(),
      shapePC1_mean = round(mean(shapePC1, na.rm=TRUE), 5),
      shapePC1_sd   = round(sd(shapePC1,   na.rm=TRUE), 5),
      shapePC2_mean = round(mean(shapePC2, na.rm=TRUE), 5),
      shapePC2_sd   = round(sd(shapePC2,   na.rm=TRUE), 5),
      across(all_of(BIO_VARS[BIO_VARS %in% names(data)]),
             list(mean=~round(mean(.x,na.rm=TRUE),4),
                  sd  =~round(sd(.x,  na.rm=TRUE),4)),
             .names="{.col}__{.fn}"),
      .groups="drop"
    ) %>%
    rename(group=1)
}

for (g in c("clade", "Order")[c("clade","Order") %in% names(sp)]) {
  tbl <- summarise_group(sp, g)
  save_csv(tbl, paste0("EFA_summary_by_", g))
}

# ── A2. PC-only summary table (matches the biomechanical 01_summary_pcs_*.csv
#        format: <group>, N, PC1_mean, PC1_sd, PC1_median, PC2_mean, PC2_sd,
#        PC2_median) — looped over every factor in GROUPS, not just clade/Order ──
summarise_pcs <- function(data, grp_col) {
  data %>%
    filter(!is.na(.data[[grp_col]])) %>%
    group_by(.data[[grp_col]]) %>%
    summarise(
      N         = n(),
      PC1_mean  = round(mean(shapePC1,   na.rm=TRUE), 3),
      PC1_sd    = round(sd(shapePC1,     na.rm=TRUE), 3),
      PC1_median= round(median(shapePC1, na.rm=TRUE), 3),
      PC2_mean  = round(mean(shapePC2,   na.rm=TRUE), 3),
      PC2_sd    = round(sd(shapePC2,     na.rm=TRUE), 3),
      PC2_median= round(median(shapePC2, na.rm=TRUE), 3),
      .groups="drop"
    ) %>%
    rename(!!grp_col := 1)
}

pcs_summary_groups <- unique(c("clade", GROUPS))
pcs_summary_groups <- pcs_summary_groups[pcs_summary_groups %in% names(sp)]

cat("\n  EFA PC summary tables (variance explained by group):\n")
for (g in pcs_summary_groups) {
  tbl <- summarise_pcs(sp, g)
  save_csv(tbl, paste0("EFA_summary_pcs_", g))
}
