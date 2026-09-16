if (exists("pruned_tree")) {
  sp$species_clean <- gsub(" ","_",trimws(sp$species))
  
  traits_phylo <- c("shapePC1","shapePC2",
                    BIO_VARS[BIO_VARS %in% names(sp)])
  
  k_results <- do.call(rbind, lapply(traits_phylo, function(tr) {
    trait_df <- sp %>% select(species_clean, value=all_of(tr)) %>% filter(!is.na(value))
    common   <- intersect(trait_df$species_clean, pruned_tree$tip.label)
    if (length(common) < 3)
      return(data.frame(trait=tr, n=length(common), K=NA, p_value=NA, sig=""))
    vals     <- setNames(trait_df$value[match(common,trait_df$species_clean)], common)
    sub_tree <- ape::drop.tip(pruned_tree, pruned_tree$tip.label[!pruned_tree$tip.label %in% common])
    vals     <- vals[sub_tree$tip.label]
    res <- tryCatch(phytools::phylosig(sub_tree, vals, method="K", test=TRUE, nsim=999),
                    error=function(e) NULL)
    if (is.null(res)) return(data.frame(trait=tr, n=length(common), K=NA, p_value=NA, sig=""))
    pv <- round(res$P,6)
    data.frame(trait=tr, n=length(common), K=round(res$K,6), p_value=pv,
               sig=case_when(pv<0.001~"***",pv<0.01~"**",pv<0.05~"*",TRUE~""))
  }))
  save_csv(k_results, "EFA_blomberg_K_all_traits")
  cat("  Blomberg's K results:\n"); print(k_results)
} else {
  cat("  pruned_tree not found — phylogenetic signal skipped.\n")
  cat("  Run the phylo section of Msc_perf_surfaces_v4.R first.\n")
}


# =============================================================================
# SECTION J — VIOLIN PLOTS: SHAPE PC BY CLADE AND SUBGROUP
# =============================================================================
cat("\n── Section J: Violin plots by clade and subgroup ──\n")

# Violin colour palette reused from CLADE_COLS if available
VIO_GROUPS <- c("clade","Diet_primary","Diet_secondary",
                "Palaeoenvironment","Depositional","Time_Bin","Diet_combo")
VIO_GROUPS <- VIO_GROUPS[VIO_GROUPS %in% names(sp)]

# Helper: stat_summary for median point + IQR vertical line (half-violin style)
add_median_stat <- function(grp_col) {
  list(
    # Vertical line spanning IQR (Q1 to Q3) at median x position
    stat_summary(aes(colour=.data[[grp_col]]),
                 fun.data=function(x) {
                   data.frame(y    = median(x, na.rm=TRUE),
                              ymin = quantile(x, 0.25, na.rm=TRUE),
                              ymax = quantile(x, 0.75, na.rm=TRUE))
                 },
                 geom="linerange", linewidth=1.4, alpha=0.90),
    # Point at median
    stat_summary(aes(colour=.data[[grp_col]]),
                 fun=function(x) median(x, na.rm=TRUE),
                 geom="point", size=3.0, shape=21,
                 fill="white", stroke=1.5)
  )
}

for (resp in c("shapePC1","shapePC2")) {
  for (grp in VIO_GROUPS) {
    d_v <- sp %>% filter(!is.na(.data[[grp]]), .data[[grp]]!="",
                         !is.na(.data[[resp]]))
    if (nrow(d_v) < 5 || length(unique(d_v[[grp]])) < 2) next
    
    ylab <- if (resp=="shapePC1") efa_lab1 else efa_lab2
    
    # KW annotation
    kw_v  <- kw_safe(d_v[[resp]], d_v[[grp]])
    
    p_vio <- ggplot(d_v, aes(fct_infreq(as.factor(.data[[grp]])),
                             .data[[resp]])) +
      # 1. Violin body
      geom_violin(aes(fill=.data[[grp]]),
                  alpha=0.45, trim=FALSE, colour="grey30",
                  linewidth=0.45) +
      # 2. Individual jitter points (small, transparent)
      geom_jitter(aes(colour=.data[[grp]]),
                  width=0.12, alpha=0.30, size=1.2) +
      # 3. Median line (IQR span) + point — half-violin style
      stat_summary(aes(colour=.data[[grp]]),
                   fun.data=function(x)
                     data.frame(y    = median(x, na.rm=TRUE),
                                ymin = quantile(x,0.25,na.rm=TRUE),
                                ymax = quantile(x,0.75,na.rm=TRUE)),
                   geom="linerange", linewidth=1.6, alpha=0.95) +
      stat_summary(aes(colour=.data[[grp]]),
                   fun=function(x) median(x, na.rm=TRUE),
                   geom="point", size=3.5, shape=21,
                   fill="white", stroke=1.8) +
      # 4. KW annotation
      annotate("text",
               x=length(unique(d_v[[grp]]))*0.5, y=Inf,
               label=if(!is.na(kw_v$H))
                 sprintf("KW: H=%.2f, p%s%.4g",
                         kw_v$H,
                         ifelse(kw_v$p_value<0.001,"<","="),
                         ifelse(kw_v$p_value<0.001,0.001,kw_v$p_value))
               else "KW: n.s.",
               hjust=0.5, vjust=1.5, size=2.8, family="mono", colour="grey30") +
      # Scales and theme
      scale_fill_brewer(palette="Pastel1", guide="none") +
      scale_colour_brewer(palette="Set1", guide="none") +
      labs(title=paste(resp, "distribution by", gsub("_"," ",grp)),
           subtitle="Violin + IQR line + median point (white) | jitter = individual specimens",
           x=NULL, y=ylab) +
      theme_bw(base_size=11) +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
            plot.title  = element_text(face="bold"),
            panel.grid.major.x = element_blank())
    
    w <- if (length(unique(d_v[[grp]])) > 10) 15 else 11
    save_plot(p_vio, paste0("EFA_violin_",resp,"_by_",grp), w=w, h=6)
  }
}