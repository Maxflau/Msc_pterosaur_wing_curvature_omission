# ── K3. MORPHOSPACE EXPANSION ~ OPTIMALITY ───────────────────────────────────
cat("\n  K3: Centroid distance (morphospace expansion) ~ mean_optimality\n")
if (nrow(time_summary) >= 4) {
  cor_exp <- cor.test(time_summary$centroid_dist, time_summary$mean_opt, method="spearman", exact=FALSE)
  lm_exp  <- lm(centroid_dist ~ mean_opt, data=time_summary)
  k3_results <- data.frame(
    Test       = "Centroid_distance ~ mean_optimality (temporal)",
    rho        = round(cor_exp$estimate, 4),
    rho_p      = round(cor_exp$p.value, 5),
    LM_R2      = round(summary(lm_exp)$r.squared, 4),
    LM_p       = round(pf(summary(lm_exp)$fstatistic[1],summary(lm_exp)$fstatistic[2],summary(lm_exp)$fstatistic[3],lower.tail=FALSE), 5)
  ) %>% mutate(sig=case_when(rho_p<0.001~"***",rho_p<0.01~"**",rho_p<0.05~"*",  rho_p<0.1~".",TRUE~""))
  save_csv(k3_results, "EFA_K3_morphospace_expansion_optimality")
  p_k3 <- ggplot(time_summary, aes(Midpoint, centroid_dist)) + geom_line(colour="#2980B9", linewidth=0.9) + geom_point(aes(fill=mean_opt, size=n_species), shape=21, colour="grey20", stroke=0.7) +
    geom_text_repel(aes(label=Time_Bin), size=2.8, colour="grey30") +
    scale_fill_viridis(name="Mean optimality", limits=c(0,1)) +
    scale_size_continuous(name="N species", range=c(3,10)) +
    scale_x_reverse(name="Age (Ma)") +
    labs(title="Morphospace expansion (centroid distance) through time",
         subtitle="Color = mean Pareto optimality per Time_Bin",
         y="Distance from overall centroid (shapePC space)") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold"))
  save_plot(p_k3, "EFA_K3_morphospace_expansion_time", w=11, h=6)
  cat(sprintf("  K3: rho=%.3f, p=%.4f %s\n",
              k3_results$rho, k3_results$rho_p, k3_results$sig))
}

# ── K4. PGLS: CLADE-LEVEL OPTIMALITY ~ CLADE RICHNESS / MORPHOLOGICAL DIVERSITY
cat("\n  K4: PGLS — clade optimality ~ species richness + morphological diversity\n")

if (exists("pruned_tree") && "mean_optimality" %in% names(sp)) {
  library(ape); library(nlme); library(phytools)
  
  # Build clade-level summary
  clade_sum <- sp %>%
    group_by(clade) %>%
    summarise(
      n_species     = n(),
      mean_opt      = mean(mean_optimality,  na.rm=TRUE),
      mean_war      = mean(aspect_ratio,     na.rm=TRUE),
      mean_svm      = mean(von_mises_stress, na.rm=TRUE),
      disparity_clade = var(shapePC1, na.rm=TRUE) + var(shapePC2, na.rm=TRUE),
      .groups       = "drop"
    ) %>%
    filter(!is.na(clade), clade != "")
  
  # Match clade names to tree tip labels (use representative species per clade)
  sp_clade_rep <- sp %>%
    group_by(clade) %>%
    slice(1) %>%
    ungroup() %>%
    select(clade, species) %>%
    mutate(tip = gsub(" ","_", trimws(species)))
  
  # Prune tree to one representative tip per clade
  keep_tips <- sp_clade_rep$tip[sp_clade_rep$tip %in% pruned_tree$tip.label]
  if (length(keep_tips) >= 4) {
    drop_tips  <- pruned_tree$tip.label[!pruned_tree$tip.label %in% keep_tips]
    clade_tree <- ape::drop.tip(pruned_tree, drop_tips)
    
    # Align clade summary to tree
    clade_df <- sp_clade_rep %>%
      filter(tip %in% clade_tree$tip.label) %>%
      left_join(clade_sum, by="clade") %>%
      filter(!is.na(mean_opt)) %>%
      as.data.frame()
    rownames(clade_df) <- clade_df$tip
    clade_df <- clade_df[clade_tree$tip.label[clade_tree$tip.label %in% rownames(clade_df)], ]
    clade_tree2 <- ape::drop.tip(clade_tree,
                                 clade_tree$tip.label[!clade_tree$tip.label %in% rownames(clade_df)])
    
    if (nrow(clade_df) >= 4) {
      # PGLS: n_species ~ mean_opt (phylogenetic generalised least squares)
      pgls_rich <- tryCatch({
        gls(n_species ~ mean_opt,
            data       = clade_df,
            correlation= corBrownian(phy=clade_tree2),
            method     = "ML")
      }, error=function(e) {cat("  PGLS rich failed:", e$message,"\n"); NULL})
      
      pgls_disp <- tryCatch({
        gls(disparity_clade ~ mean_opt,
            data       = clade_df,
            correlation= corBrownian(phy=clade_tree2),
            method     = "ML")
      }, error=function(e) {cat("  PGLS disp failed:", e$message,"\n"); NULL})
      
      k4_results <- data.frame()
      for (obj in list(pgls_rich, pgls_disp)) {
        if (!is.null(obj)) {
          s    <- summary(obj)
          resp <- as.character(formula(obj)[[2]])
          coef <- s$tTable["mean_opt",]
          k4_results <- rbind(k4_results, data.frame(
            Response  = resp,
            Predictor = "mean_optimality",
            Estimate  = round(coef[1],6),
            SE        = round(coef[2],6),
            t_value   = round(coef[3],3),
            p_value   = round(coef[4],5),
            n_clades  = nrow(clade_df),
            Method    = "PGLS (Brownian motion)"
          ))
        }
      }
      if (nrow(k4_results) > 0) {
        k4_results <- k4_results %>%
          mutate(sig=case_when(p_value<0.001~"***",p_value<0.01~"**",
                               p_value<0.05~"*",  p_value<0.1~".",TRUE~""))
        save_csv(k4_results, "EFA_K4_PGLS_optimality_richness_diversity")
        cat("  K4 PGLS results:\n"); print(k4_results)
      }
    }
  } else {
    cat("  <4 tree tips matched — PGLS skipped.\n")
  }
} else {
  cat("  pruned_tree or mean_optimality not available — K4 PGLS skipped.\n")
}