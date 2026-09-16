PTERODACT <- c("Azhdarchoidea","Ctenochasmatoidea","Dsungaripteroidea",
               "Ornithocheiromorpha","Pteranodontia","basal pterodactyloidea")

if (!file.exists("output/phylogenetics/Henry_updated.nex")) {
  stop("Henry_updated.nex not found - phylomorphospace cannot proceed.")
}

all_trees  <- ape::read.nexus("output/phylogenetics/Henry_updated.nex")
ptero_tree <- if (inherits(all_trees, "multiPhylo")) all_trees[[1]] else all_trees
cat(sprintf("Tree: %d tips, %d nodes\n", ape::Ntip(ptero_tree), ape::Nnode(ptero_tree)))
tip_clean <- gsub("_", " ", ptero_tree$tip.label)
match_lookup <- data.frame(species=character(), tip=character(), stringsAsFactors=FALSE)
for (sp in shape_perf$species) {
  sp_t <- trimws(sp)
  ex   <- which(tip_clean == sp_t)
  if (length(ex) > 0) {
    match_lookup <- rbind(match_lookup, data.frame(species=sp, tip=ptero_tree$tip.label[ex[1]]))
    next
  }
  parts <- strsplit(sp_t, " ")[[1]]
  if (length(parts) >= 2) {
    pm <- which(grepl(paste(parts[1], parts[2]), tip_clean, ignore.case=TRUE))
    if (length(pm) > 0)
      match_lookup <- rbind(match_lookup, data.frame(species=sp, tip=ptero_tree$tip.label[pm[1]]))
  }
}
cat(sprintf("Matched %d / %d species to tree\n", nrow(match_lookup), nrow(shape_perf)))

phylo_sp <- shape_perf %>%
  inner_join(match_lookup, by="species") %>%
  filter(!duplicated(tip), !is.na(shapePC1), !is.na(shapePC2), !is.na(pareto_rank_ratio))

tips_keep <- phylo_sp$tip[phylo_sp$tip %in% ptero_tree$tip.label]
pruned    <- ape::drop.tip(ptero_tree, ptero_tree$tip.label[!ptero_tree$tip.label %in% tips_keep])
cat(sprintf("Pruned tree: %d tips\n", ape::Ntip(pruned)))

spc1 <- setNames(phylo_sp$shapePC1[match(pruned$tip.label, phylo_sp$tip)], pruned$tip.label)
spc2 <- setNames(phylo_sp$shapePC2[match(pruned$tip.label, phylo_sp$tip)], pruned$tip.label)
opt  <- setNames(phylo_sp$pareto_rank_ratio[match(pruned$tip.label, phylo_sp$tip)], pruned$tip.label)

cat("Estimating ancestral states...\n")
anc_pc1 <- phytools::fastAnc(pruned, spc1)
anc_pc2 <- phytools::fastAnc(pruned, spc2)
all_pc1 <- c(spc1, anc_pc1); all_pc2 <- c(spc2, anc_pc2)
all_opt <- c(opt, phytools::fastAnc(pruned, opt))

n_tips <- ape::Ntip(pruned)
edge_df <- lapply(seq_len(nrow(pruned$edge)), function(i) {
  par <- pruned$edge[i,1]; ch <- pruned$edge[i,2]
  p_nm <- as.character(par)
  c_nm <- if (ch <= n_tips) pruned$tip.label[ch] else as.character(ch)
  if (!p_nm %in% names(all_pc1) || !c_nm %in% names(all_pc1)) return(NULL)
  data.frame(x1=all_pc1[p_nm], y1=all_pc2[p_nm], x2=all_pc1[c_nm], y2=all_pc2[c_nm])
}) %>% bind_rows()
cat(sprintf("Created %d phylogenetic edges\n", nrow(edge_df)))

K1 <- phytools::phylosig(pruned, spc1, method="K", test=TRUE, nsim=999)
K2 <- phytools::phylosig(pruned, spc2, method="K", test=TRUE, nsim=999)
cat(sprintf("K(shapePC1)=%.3f p=%.3f | K(shapePC2)=%.3f p=%.3f\n", K1$K, K1$P, K2$K, K2$P))

k_annot <- sprintf(
  "Blomberg's K\nshapePC1: K=%.3f  p%s%.3f\nshapePC2: K=%.3f  p%s%.3f",
  K1$K, ifelse(K1$P<0.001,"<","="), max(K1$P, 0.001),
  K2$K, ifelse(K2$P<0.001,"<","="), max(K2$P, 0.001)
)

make_phylo_v4 <- function(subset=c("all","non_pterodact","pterodact",
                                   "non_pterodactyliformes")) {
  subset <- match.arg(subset)
  
  # Pterodactyliformes = Darwinoptera + Pterodactyloidea (sensu Unwin 2003)
  PTERODACTYLIFORMES <- c(PTERODACT, "Darwinoptera")
  
  all_pts <- phylo_sp %>%
    mutate(grade = case_when(
      clade %in% PTERODACT        ~ "Pterodactyloidea",
      clade == "Darwinoptera"     ~ "Darwinoptera (Pterodactyliform)",
      TRUE                        ~ "Non-Pterodactyliformes"
    ))
  
  fg <- switch(subset,
               all                    = all_pts,
               non_pterodact          = filter(all_pts, !clade %in% PTERODACT),
               pterodact              = filter(all_pts,  clade %in% PTERODACT),
               non_pterodactyliformes = filter(all_pts, !clade %in% PTERODACTYLIFORMES)
  )
  
  ttl <- switch(subset,
                all                    = "Pterosauria Phylomorphospace",
                non_pterodact          = "(a) Non-Pterodactyloidea Phylomorphospace",
                pterodact              = "(b) Pterodactyloidea Phylomorphospace",
                non_pterodactyliformes = "(a) Non-Pterodactyliformes Phylomorphospace"
  )
  sub <- switch(subset,
                all                    = "All taxa — Pareto optimality landscape",
                non_pterodact          = "Basal pterosaurs projected onto optimality landscape",
                pterodact              = "Derived pterosaurs projected onto optimality landscape",
                non_pterodactyliformes = "Anurognathidae, Rhamphorhynchidae, Non-breviquartossan,\nBasal pterosaur — Pareto optimality landscape"
  )
  
  xrng <- range(all_pts$shapePC1, na.rm=TRUE)
  yrng <- range(all_pts$shapePC2, na.rm=TRUE)
  
  p <- ggplot()
  if (!is.null(pareto_bg))
    p <- p +
    geom_raster(data=pareto_bg, aes(x,y,fill=z), interpolate=TRUE) +
    scale_fill_gradient(low="white", high="#1A5C2A", name="Pareto\nOptimality",
                        limits=c(0,1), breaks=c(0,.25,.5,.75,1),
                        guide=guide_colorbar(order=99, barheight=unit(3.5,"cm"), barwidth=unit(0.32,"cm"), title.theme=element_text(size=7, face="bold")))
  if (nrow(impossible_region) > 0)
    p <- p + geom_tile(data=impossible_region, aes(x,y),
                       fill="grey80", alpha=0.75, inherit.aes=FALSE)
  p <- p +
    new_scale_fill() +
    geom_polygon(data=grid_dense, aes(gx,gy,group=grid_id),
                 fill="grey88", colour="grey55", linewidth=0.06, alpha=0.38) +
    geom_segment(data=edge_df, aes(x=x1,y=y1,xend=x2,yend=y2),
                 colour="grey30", linewidth=0.30, alpha=0.55)
  if (subset != "all")
    p <- p + geom_point(data=all_pts, aes(shapePC1,shapePC2),
                        shape=21, fill="grey75", colour="grey50",
                        size=1.5, stroke=0.25, alpha=0.30)
  
  p <- p +
    new_scale_fill() +
    geom_point(data=fg, aes(shapePC1, shapePC2),
               shape=21, fill="black", colour="black",
               size=2.4, stroke=0.4) +
    annotate("text", x=xrng[1]+diff(xrng)*0.02, y=yrng[2]-diff(yrng)*0.03,
             label=k_annot, hjust=0, vjust=1, size=2.6, colour="grey15", family="mono") +
    labs(title=ttl, subtitle=sub, x=lab1, y=lab2) +
    coord_equal(expand=FALSE) + base_theme
  p
}