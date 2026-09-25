# Test whether each performance/shape variable carries a "phylogenetic signal"
# - i.e. whether closely related species tend to have similar values, which a
# plain (non-phylogenetic) statistical test would otherwise ignore.
#
# Two measures are reported side by side:
#   - Pagel's lambda (primary measure): 0 = no signal, 1 = signal matches
#     exactly what the tree predicts under simple Brownian motion.
#   - Blomberg's K (robustness check): < 1 = less signal than the tree
#     predicts, > 1 = more signal than the tree predicts.

cat("\n================================================================================\n")
cat("PHYLOGENETIC SIGNAL - Pagel's lambda and Blomberg's K\n")
cat("================================================================================\n\n")

# --------------------------------------------------------------------------------
# 1. Find the tree object automatically, by its CLASS ("phylo"), instead of by
# a hard-coded name - so this script keeps working even if the tree variable
# was renamed earlier in the pipeline.
# --------------------------------------------------------------------------------
env_objects <- ls(envir = .GlobalEnv)
is_phylo    <- sapply(env_objects, function(nm) inherits(get(nm, envir = .GlobalEnv), "phylo"))
phylo_objects <- env_objects[is_phylo]

if (length(phylo_objects) == 0) {
  stop("No object of class 'phylo' found in the global environment.\n",
       "  Source the phylomorphospace script first.")
}

# If several tree objects exist, pick the smallest one (fewest tips) - this is
# almost always the tree that has already been pruned to match the data.
n_tips <- sapply(phylo_objects, function(nm) length(get(nm, envir = .GlobalEnv)$tip.label))
if (length(phylo_objects) > 1) {
  cat("Multiple 'phylo' objects found:\n")
  print(data.frame(object = phylo_objects, n_tips = n_tips))
  tree_name <- phylo_objects[which.min(n_tips)]
  cat(sprintf("\nUsing '%s' (fewest tips). Override tree_name manually if wrong.\n\n", tree_name))
} else {
  tree_name <- phylo_objects[1]
}

tree <- get(tree_name, envir = .GlobalEnv)
cat(sprintf("Using tree object '%s': %d tips.\n", tree_name, length(tree$tip.label)))

# --------------------------------------------------------------------------------
# 2. Match species names (performance_data_clean) to the tree's tip labels.
# Small text differences (extra spaces, underscores vs spaces) are the most
# common reason a species fails to match, so names are normalised first.
# --------------------------------------------------------------------------------
norm_name <- function(x) gsub("[[:space:]*]+", "*", trimws(x))

tip_norm <- norm_name(tree$tip.label)                       # same order/length as tree$tip.label
sp_norm  <- norm_name(performance_data_clean$species)

match_idx <- match(sp_norm, tip_norm)   # position in tree$tip.label for each species, or NA
n_matched <- sum(!is.na(match_idx))
cat(sprintf("Species matched to tree tips (after name normalisation): %d / %d\n",
            n_matched, length(sp_norm)))

if (n_matched < 4) {
  cat("\nExample tree tip labels:\n"); print(head(tree$tip.label, 10))
  cat("\nExample unmatched species names:\n")
  print(head(performance_data_clean$species[is.na(match_idx)], 10))
  stop("Fewer than 4 species matched - check naming convention above.")
}

# Direct species -> original tip label lookup (no double indirection).
tip_label_for_species <- rep(NA_character_, length(sp_norm))
tip_label_for_species[!is.na(match_idx)] <- tree$tip.label[match_idx[!is.na(match_idx)]]
species_to_tip <- setNames(tip_label_for_species, performance_data_clean$species)
species_to_tip <- species_to_tip[!is.na(species_to_tip)]

# --------------------------------------------------------------------------------
# 3. Run the signal test for every variable of interest.
# --------------------------------------------------------------------------------
SIGNAL_VARS <- intersect(
  c("aspect_ratio", "r2_hat", "von_mises_stress", "wing_curvature",
    "shape_complexity", "wing_loading_ratio", "PC1", "PC2"),
  colnames(performance_data_clean)
)

phylo_signal_results <- lapply(SIGNAL_VARS, function(v) {

  dat <- performance_data_clean[[v]]
  names(dat) <- performance_data_clean$species
  dat <- dat[is.finite(dat)]

  usable_species <- intersect(names(dat), names(species_to_tip))
  if (length(usable_species) < 4) {
    return(data.frame(variable = v, n = length(usable_species),
                      lambda = NA, lambda_p = NA, K = NA, K_p = NA))
  }

  dat_v <- dat[usable_species]
  names(dat_v) <- species_to_tip[usable_species]   # switch names to the tree's own tip labels

  tree_v <- keep.tip(tree, unname(names(dat_v)))
  if (is.null(tree_v) || length(tree_v$tip.label) < 4) {
    return(data.frame(variable = v, n = length(usable_species),
                      lambda = NA, lambda_p = NA, K = NA, K_p = NA))
  }

  sig_lambda <- phylosig(tree_v, dat_v, method = "lambda", test = TRUE)
  sig_K      <- phylosig(tree_v, dat_v, method = "K",      test = TRUE, nsim = 999)

  data.frame(
    variable = v,
    n        = length(dat_v),
    lambda   = unname(sig_lambda$lambda),
    lambda_p = unname(sig_lambda$P),
    K        = unname(sig_K$K),
    K_p      = unname(sig_K$P)
  )
})

phylo_signal_results <- dplyr::bind_rows(phylo_signal_results)

cat("\nPhylogenetic signal - Pagel's lambda (primary) vs Blomberg's K (robustness check):\n")
phylo_signal_display <- phylo_signal_results
phylo_signal_display$lambda   <- format_fixed(phylo_signal_results$lambda, 3)
phylo_signal_display$lambda_p <- format_fixed(phylo_signal_results$lambda_p, 3)
phylo_signal_display$K        <- format_fixed(phylo_signal_results$K, 3)
phylo_signal_display$K_p      <- format_fixed(phylo_signal_results$K_p, 3)
print(phylo_signal_display, row.names = FALSE)

# Flag any variable where the two measures disagree on significance (p < 0.05) -
# worth a closer look before trusting either one on its own.
phylo_signal_results$disagreement <-
  (phylo_signal_results$lambda_p < 0.05) != (phylo_signal_results$K_p < 0.05)

if (any(phylo_signal_results$disagreement, na.rm = TRUE)) {
  cat("\nVariables where lambda and K disagree on significance (p < 0.05):\n")
  print(phylo_signal_results[which(phylo_signal_results$disagreement),
                             c("variable", "lambda", "lambda_p", "K", "K_p")])
}

dir.create("output/results/supplementals", showWarnings = FALSE, recursive = TRUE)
write.csv(phylo_signal_results, "output/results/supplementals/phylogenetic_signal_lambda_K.csv",
          row.names = FALSE)

cat("\nAvailable: phylo_signal_results, tree, species_to_tip\n\n")
