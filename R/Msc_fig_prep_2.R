# Set clade colours and build a smooth background stress surface.
# Use one fixed border colour for each clade.
CLADE_COLS <- c(
  "Anurognathidae"         = "#46852F",
  "Azhdarchoidea"          = "#DB8E00",
  "Basal pterodactyliform" = "#F4A7B9",
  "Basal pterosaur"        = "#0CE27F",
  "Ctenochasmatoidea"      = "#9F1212",
  "Darwinoptera"           = "#791444",
  "Dsungaripteroidea"      = "#60A293",
  "Eudimorphodontoidea"    = "#FF00F0",
  "Ornithocheiromorpha"    = "#FFE413",
  "Pteranodontia"          = "#DA8AAA",
  "Rhamphorhynchidae"      = "#03240B"
)

match_clade_colour <- function(clade_names) {
  norm <- function(x) tolower(trimws(x))
  lut <- CLADE_COLS
  names(lut) <- norm(names(CLADE_COLS))
  out <- unname(lut[norm(clade_names)])
  names(out) <- clade_names
  out
}

# Make sure the figure folders exist before saving plots.
for (d in c("output/plots/biomechanic_plot",
            "output/plots_PDF/biomechanic_plot")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("stress_grid")) stop("stress_grid not found - source Msc_plot_layer_v3.R first.")
if (!exists("wing_bg"))     stop("wing_bg not found - source Msc_wing_overlay_v1.R first.")

# Build the axis label text from the current PCA output.
pc_lab <- function(k) {
  sprintf("PC%d (%.1f%%)", k, summary(pca_performance)$importance[2, k] * 100)
}
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required for the smoothed stress surface: install.packages(\"mgcv\")")
}

# Smooth the stress values across the morphospace grid for the background layer.
build_smooth_stress_grid <- function(value_col = "von_mises_stress") {
  df <- data.frame(x = performance_data_clean$PC1, y = performance_data_clean$PC2,
                   z = performance_data_clean[[value_col]])
  df <- df[complete.cases(df) & is.finite(df$z) & df$z > 0, ]
  if (nrow(df) < 15) { cat("  Too few points for a smoothed stress surface.\n"); return(NULL) }
  m <- mgcv::gam(log(z) ~ s(x, y, k = min(60, nrow(df) - 1)), data = df)
  g <- data.frame(PC1 = stress_grid$PC1, PC2 = stress_grid$PC2)
  g$stress <- exp(predict(m, newdata = data.frame(x = g$PC1, y = g$PC2)))
  g
}

stress_grid_smooth <- build_smooth_stress_grid("von_mises_stress")
if (!is.null(stress_grid_smooth)) {
  cat(sprintf("Smoothed stress surface (GAM) built: %d cells, range %.4g to %.4g\n\n",
              nrow(stress_grid_smooth),
              min(stress_grid_smooth$stress), max(stress_grid_smooth$stress)))
}
