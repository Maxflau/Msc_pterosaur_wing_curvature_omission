cat("\n================================================================================\n")
cat("Plot preparations\n")
cat("================================================================================\n\n")

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  stop("Package 'ggnewscale' is required: install.packages(\"ggnewscale\")")
}

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

# Used for the "Environment" column
ENV_COLS <- c(
  "Archipelago"                   = "#A300E9",
  "Coastal environment"           = "#4F4F4F",
  "Desert"                        = "#B51963",
  "Floodplain"                    = "#FFC107",
  "Fluvial environment"           = "#FF5722",
  "Forested wetland environment"  = "#FF9800",
  "Marine environment"            = "#FAAF90",
  "Semi-arid floodplain"          = "#4CAF50"
)

DIET_COLS <- c(
  "Carnivore"      = "#E74C3C",
  "Durophageous"   = "#4CAF50",
  "Filter Feeder"  = "#7B1FA2",
  "Generalist"     = "#144B00",
  "Herbivore"      = "#80DEEA",
  "Insectivore"    = "#FFA726",
  "Piscivore"      = "#F5F5B8"
)

DIETCOMBO_COLS <- c(
  "Carnivore+Durophageous"   = "#9575CD", "Carnivore+Generalist"     = "#B0BEC5",
  "Carnivore+Insectivore"    = "#42A5F5", "Carnivore+Piscivore"      = "#7B1FA2",
  "Durophageous+Piscivore"   = "#FF8A65", "Filter Feeder+Piscivore"  = "#FF7043",
  "Generalist+Piscivore"     = "#EC407A", "Herbivore+Durophageous"   = "#F48FB1",
  "Insectivore+Carnivore"    = "#26A69A", "Insectivore+Durophageous" = "#CE93D8",
  "Insectivore+Piscivore"    = "#EF5350", "Piscivore+Carnivore"      = "#FFA726",
  "Piscivore+Durophageous"   = "#66BB6A", "Piscivore+Generalist"     = "#FFEE58",
  "Piscivore+Insectivore"    = "#FFC107", "Piscivore+Piscivore"      = "#FDD835"
)

# Used for the "Depositional.settings.paleoenvironment" column
DEPOSITIONAL_COLS <- c(
  "Aeolian"                 = "#482072",
  "Alluvial plain"          = "#4CAF50",
  "Coastal/Shallow marine"  = "#D55E00",
  "Fluviodeltaic"           = "#144B00",
  "Lacustrine: large lake"  = "#FFC107",
  "Lacustrine: small lake"  = "#F0CD5D",
  "Lagoonal deposit"        = "#A300E9",
  "playa"                   = "#6A9DA3"
)

# Generic case/whitespace-insensitive lookup against any fixed colour set,
# so labelling differences in the data (capitalisation, extra spaces around
# separators like " + ") don't silently fall through.
match_fixed_colour <- function(names_in_data, colour_lut) {
  norm <- function(x) gsub("\\s+", "", tolower(trimws(x)))
  lut <- colour_lut
  names(lut) <- norm(names(colour_lut))
  out <- unname(lut[norm(names_in_data)])
  names(out) <- names_in_data
  out
}

for (d in c("clade_analyses_v4/plots/biomechanic_plot",
            "clade_analyses_v4/plots_PDF/biomechanic_plot")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("stress_grid")) stop("stress_grid not found - source Msc_plot_layer_v3.R first.")
if (!exists("wing_bg"))     stop("wing_bg not found - source Msc_wing_overlay_v1.R first.")

pc_lab <- function(k) {
  sprintf("PC%d (%.1f%%)", k, summary(pca_performance)$importance[2, k] * 100)
}

METRIC_LABELS <- c(
  von_mises_stress = "Von Mises stress\n(scale-invariant)",r2_hat = "Second moment\nof area", aspect_ratio = "Aspect ratio", wing_curvature = "Wing curvature",
  shape_complexity = "Shape complexity", pareto_rank_ratio = "Pareto\nrank ratio")

surf_label <- if (BACKGROUND_VAR %in% names(METRIC_LABELS)) {
  METRIC_LABELS[[BACKGROUND_VAR]]
} else BACKGROUND_VAR
cat("Surface variable:", BACKGROUND_VAR, "| colour bar reads:",gsub("\n", " ", surf_label), "\n\n")
