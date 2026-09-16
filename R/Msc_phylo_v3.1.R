cat("\n================================================================================\n")
cat("PHYLOGENETIC PREPARATION\n")
cat("================================================================================\n\n")

# Load the tree packages used in the next steps.
for (pkg in c("ape", "phytools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.")
  }
}
library(ape)
library(phytools)

TREE_FILE <- "output/phylogenetics/Henry_updated.nex"

##################################################################################
# 1. Load the tree used for the phylogenetic analyses.
##################################################################################

if (!file.exists(TREE_FILE)) {
  stop("Tree file '", TREE_FILE, "' not found.")
}

all_trees <- ape::read.nexus(TREE_FILE)

if (inherits(all_trees, "multiPhylo")) {
  ptero_tree <- all_trees[[1]]
  n_trees <- length(all_trees)
  cat(sprintf("Loaded %d trees; using tree 1 for the main analysis.\n", n_trees))
} else {
  ptero_tree <- all_trees
  all_trees <- NULL
  n_trees <- 1
  cat("Loaded a single tree.\n")
}

if (is.null(ptero_tree$edge.length)) {
  stop("Tree has no branch lengths — PGLS and signal tests require them.")
}
cat(sprintf("Tree: %d tips, %d internal nodes\n\n",Ntip(ptero_tree), Nnode(ptero_tree)))

# ####################################################################################
# 2. Match specimen names to the tree tip names.
# ####################################################################################

# Clean names in the same way on both sides before matching them.
normalise_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[ ]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  tolower(x)
}

phylo_data <- performance_data_clean
phylo_data$tip_key <- normalise_name(phylo_data$species)
tip_key <- normalise_name(ptero_tree$tip.label)

matched <- phylo_data$tip_key %in% tip_key
cat(sprintf("Specimens matched to tips: %d / %d\n", sum(matched), nrow(phylo_data)))

if (any(!matched)) {
  cat("\nUnmatched specimens (check for synonymy or spelling before accepting):\n")
  print(phylo_data$species[!matched])
  cat("\n")
}
if (sum(matched) < 20) {
  cat("First 10 tip labels for comparison:\n")
  print(head(ptero_tree$tip.label, 10))
  stop("Too few matches — resolve naming before proceeding.")
}

phylo_data <- phylo_data[matched, ]
phylo_data <- phylo_data[!duplicated(phylo_data$tip_key), ]

# ####################################################################################
# 3. Drop unused tips, then line the data up with the pruned tree.
# ####################################################################################
keep_tips <- ptero_tree$tip.label[tip_key %in% phylo_data$tip_key]
phy_pruned <- drop.tip(ptero_tree, setdiff(ptero_tree$tip.label, keep_tips))
rownames(phylo_data) <- phylo_data$tip_key
phylo_data <- phylo_data[normalise_name(phy_pruned$tip.label), ]
stopifnot(all(normalise_name(phy_pruned$tip.label) == phylo_data$tip_key))

# Zero-length branches break the variance-covariance matrix
if (any(phy_pruned$edge.length <= 0)) {
  n_zero <- sum(phy_pruned$edge.length <= 0)
  phy_pruned$edge.length[phy_pruned$edge.length <= 0] <-
    1e-6 * max(phy_pruned$edge.length)
  cat(sprintf("Adjusted %d zero-length branches.\n", n_zero))
}

cat(sprintf("Pruned tree: %d tips; analysis dataset: %d rows\n\n",
            Ntip(phy_pruned), nrow(phylo_data)))

# Order group, used by every downstream plot. The database's Order field is
# trusted for everything EXCEPT Darwinoptera: it sits outside Pterodactyloidea
# (Unwin 2003) but is not consistently labelled "Non Pterodactyliform" in the
# raw data, so it is forced into that side explicitly rather than left to
# whatever the database happens to say.
phylo_data$order_group <- ifelse(
  phylo_data$Order == "Non Pterodactyliform" | phylo_data$clade == "Darwinoptera",
  "Non-pterodactyliform", "pterodactyloidea")

cat("Specimens per order group:\n")
print(table(phylo_data$order_group))
cat("\n")

#####################################################################################
# 4. Ancestral states and phylomorphospace edges
#####################################################################################
pc1_vals <- setNames(phylo_data$PC1, phy_pruned$tip.label)
pc2_vals <- setNames(phylo_data$PC2, phy_pruned$tip.label)

cat("Estimating ancestral states...\n")
anc_pc1 <- phytools::fastAnc(phy_pruned, pc1_vals)
anc_pc2 <- phytools::fastAnc(phy_pruned, pc2_vals)

all_pc1 <- c(pc1_vals, anc_pc1)
all_pc2 <- c(pc2_vals, anc_pc2)

# Build the branch table used later for tree overlays.
n_tips <- Ntip(phy_pruned)
parent <- phy_pruned$edge[, 1]
child  <- phy_pruned$edge[, 2]
child_name  <- ifelse(child <= n_tips, phy_pruned$tip.label[child], as.character(child))
parent_name <- as.character(parent)
edge_data <- data.frame(
  x_start = all_pc1[parent_name], y_start = all_pc2[parent_name],
  x_end   = all_pc1[child_name],  y_end   = all_pc2[child_name])
edge_data <- edge_data[complete.cases(edge_data), ]

cat(sprintf("Phylogenetic edges: %d\n\n", nrow(edge_data)))

cat("Available: phy_pruned, phylo_data, pc1_vals, pc2_vals, edge_data, all_trees\n\n")
