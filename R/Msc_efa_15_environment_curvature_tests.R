# Run extra group tests, then check environment and curvature effects.
kw_safe <- function(vals, grps) {
  ok <- !is.na(vals) & !is.na(grps) & grps != ""
  vals <- vals[ok]; grps <- droplevels(as.factor(as.character(grps[ok])))
  if (nlevels(grps) < 2 || length(vals) < 4) return(data.frame(H=NA,df_kw=NA,p_value=NA,sig=""))
  res <- tryCatch(kruskal.test(vals, grps), error=function(e) NULL)
  if (is.null(res)) return(data.frame(H=NA,df_kw=NA,p_value=NA,sig=""))
  pv <- round(res$p.value,6)
  data.frame(H=round(res$statistic,3), df_kw=res$parameter, p_value=pv,
             sig=case_when(pv<0.001~"***",pv<0.01~"**",pv<0.05~"*",TRUE~""))
}

kw_results <- do.call(rbind, lapply(avail_groups, function(g)
  do.call(rbind, lapply(c("shapePC1","shapePC2"), function(v) {
    if (!g %in% names(sp)) return(NULL)
    res <- kw_safe(sp[[v]], sp[[g]])
    cbind(Group=g, Variable=v, res)
  }))
))
save_csv(kw_results, "EFA_kruskal_wallis_shapePC")

# =============================================================================
# SECTION F — ENVIRONMENTAL CONTROL (Msc_code_environment_control_v2)
# KW for shapePC2 vs preservation quality
# Chi-square: habitat vs depositional setting
# =============================================================================
# Test whether preservation and setting help explain shape differences.
cat("\n── Section F: Environmental control ──\n")

if ("Depositional" %in% names(sp) && !all(is.na(sp$Depositional))) {
  env_raw <- trimws(sp$Depositional)
  sp$Pres_Quality <- case_when(
    grepl("Lagoon",                 env_raw, ignore.case=TRUE) ~ "High",
    grepl("Lacustrine.*small",      env_raw, ignore.case=TRUE) ~ "High",
    grepl("Lacustrine.*large",      env_raw, ignore.case=TRUE) ~ "Med-High",
    grepl("Fluviodeltaic",          env_raw, ignore.case=TRUE) ~ "Medium",
    grepl("Coastal|Shallow marine", env_raw, ignore.case=TRUE) ~ "Medium",
    grepl("Alluvial",               env_raw, ignore.case=TRUE) ~ "Low",
    grepl("playa",                  env_raw, ignore.case=TRUE) ~ "Low",
    grepl("Aeolian",                env_raw, ignore.case=TRUE) ~ "Low",
    TRUE ~ NA_character_
  )
  sp$Pres_Quality <- factor(sp$Pres_Quality, levels=c("Low","Medium","Med-High","High"))

  sp_env <- sp %>% filter(!is.na(Pres_Quality))
  kw_pq <- kruskal.test(shapePC2 ~ Pres_Quality, data=sp_env)
  pq_summary <- sp_env %>%
    group_by(Pres_Quality) %>%
    summarise(n=n(), mean_PC2=round(mean(shapePC2,na.rm=TRUE),5),
              sd_PC2=round(sd(shapePC2,na.rm=TRUE),5), .groups="drop")
  pq_summary$KW_H  <- round(kw_pq$statistic, 3)
  pq_summary$KW_p  <- round(kw_pq$p.value, 5)
  save_csv(pq_summary, "EFA_preservation_quality_shapePC2")

  p_pq <- ggplot(sp_env, aes(Pres_Quality, shapePC2, fill=Pres_Quality)) +
    geom_boxplot(alpha=0.55, outlier.shape=NA) +
    geom_jitter(width=0.15, alpha=0.4, size=1.5) +
    scale_fill_brewer(palette="RdYlGn") +
    annotate("text", x=2.5, y=max(sp_env$shapePC2,na.rm=TRUE)*0.95,
             label=sprintf("KW: H=%.2f, p=%.4f", kw_pq$statistic, kw_pq$p.value),
             size=3.5, family="mono") +
    labs(title="EFA shapePC2 by Preservation Quality",
         subtitle="Taphonomic bias control (adapted from Msc_depo_stats.R)",
         x="Preservation quality", y=efa_lab2) +
    theme_bw() + theme(legend.position="none")
  save_plot(p_pq, "EFA_preservation_quality_shapePC2")
}

# Chi-square: habitat × depositional setting
if ("Palaeoenvironment" %in% names(sp) && "Depositional" %in% names(sp)) {
  sp_ct <- sp %>% filter(!is.na(Palaeoenvironment), Palaeoenvironment!="",
                         !is.na(Depositional),       Depositional!="")
  ct  <- table(sp_ct$Palaeoenvironment, sp_ct$Depositional)
  chi <- chisq.test(ct, simulate.p.value=TRUE, B=10000)
  cv  <- sqrt(as.numeric(chi$statistic) / (sum(ct)*(min(dim(ct))-1)))
  chi_tbl <- data.frame(X2=round(chi$statistic,3), p_value=round(chi$p.value,4),
                        CramersV=round(cv,3), B_sims=10000)
  save_csv(chi_tbl, "EFA_env_depo_chisquare")
  cat(sprintf("  Chi-square: X2=%.2f, p=%.4f, Cramer's V=%.3f\n",
              chi$statistic, chi$p.value, cv))
}

# =============================================================================
# SECTION G — CURVATURE DISTRIBUTION ANALYSIS (Msc_curvature_distri.R)
# Wing curvature presence vs shapePC1/PC2 positions
# =============================================================================
# Compare specimens with and without measured wing curvature.
cat("\n── Section G: Curvature distribution ──\n")

if ("wing_curvature" %in% names(sp)) {
  sp$has_curvature <- !is.na(sp$wing_curvature) & sp$wing_curvature > 0

  curv_summary <- sp %>%
    group_by(has_curvature) %>%
    summarise(n=n(),
              mean_shapePC1 = round(mean(shapePC1,na.rm=TRUE),5),
              sd_shapePC1   = round(sd(shapePC1,  na.rm=TRUE),5),
              mean_shapePC2 = round(mean(shapePC2,na.rm=TRUE),5),
              sd_shapePC2   = round(sd(shapePC2,  na.rm=TRUE),5),
              mean_WAR      = round(mean(aspect_ratio,na.rm=TRUE),3),
              .groups="drop")
  save_csv(curv_summary, "EFA_curvature_presence_summary")

  p_c1 <- ggplot(sp, aes(has_curvature, shapePC1, fill=has_curvature)) +
    geom_boxplot(alpha=0.55, outlier.shape=NA) + geom_jitter(width=0.15,alpha=0.4,size=1.5) +
    scale_fill_viridis(discrete=TRUE, guide="none") +
    labs(title="EFA shapePC1 by wing curvature presence", x="Curvature measured?", y=efa_lab1) +
    theme_bw()
  p_c2 <- ggplot(sp, aes(has_curvature, shapePC2, fill=has_curvature)) +
    geom_boxplot(alpha=0.55, outlier.shape=NA) + geom_jitter(width=0.15,alpha=0.4,size=1.5) +
    scale_fill_viridis(discrete=TRUE, guide="none") +
    labs(title="EFA shapePC2 by wing curvature presence", x="Curvature measured?", y=efa_lab2) +
    theme_bw()
  p_c3 <- ggplot(sp, aes(has_curvature, aspect_ratio, fill=has_curvature)) +
    geom_boxplot(alpha=0.55, outlier.shape=NA) + geom_jitter(width=0.15,alpha=0.4,size=1.5) +
    scale_fill_viridis(discrete=TRUE, guide="none") +
    labs(title="WAR by curvature presence", x="Curvature measured?", y="Aspect Ratio") +
    theme_bw()
  p_curv <- ggarrange(p_c1, p_c2, p_c3, ncol=3, labels=c("A","B","C"))
  save_plot(p_curv, "EFA_curvature_presence_sensitivity", w=15, h=5)
}
