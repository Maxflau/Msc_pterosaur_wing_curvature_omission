if (!exists("phylo_sp") || !exists("make_phylo_v4")) {
  stop("phylo_sp / make_phylo_v4() not found - source Msc_EFA_2_6a.R first.")
}

violin_df <- phylo_sp %>%
  mutate(grade = case_when(
    clade %in% PTERODACT    ~ "Pterodactyliformes",
    clade == "Darwinoptera" ~ "Pterodactyliformes",
    TRUE                    ~ "Non-Pterodactyliformes"
  ))
wx <- wilcox.test(pareto_rank_ratio ~ grade, data=violin_df)
cat(sprintf("Wilcoxon (Non-Pterodactyliformes vs Pterodactyliformes): W=%.0f  p=%.3e\n",
            wx$statistic, wx$p.value))

p_violin <- ggplot(violin_df, aes(grade, pareto_rank_ratio, fill=grade)) +
  geom_violin(alpha=0.55, trim=FALSE, colour="grey30") +
  geom_boxplot(width=0.12, fill="white", alpha=0.85, outlier.shape=NA) +
  geom_jitter(width=0.08, alpha=0.30, size=1.5, colour="grey30") +
  scale_fill_manual(values=c("Non-Pterodactyliformes"="#4292c6",
                             "Pterodactyliformes"="#ef6548"), guide="none") +
  annotate("text", x=1.5, y=max(violin_df$pareto_rank_ratio, na.rm=TRUE)*0.97,
           label=sprintf("Wilcoxon\nW=%.0f\np=%.2e", wx$statistic, wx$p.value),
           size=3, hjust=0.5, family="mono") +
  labs(title="(c) Pareto Optimality by Grade",
       subtitle="Non-Pterodactyliformes vs Pterodactyliformes",
       x="", y="Pareto Optimality") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))

traits_K <- c("aspect_ratio","r2_hat","wing_loading","von_mises_stress","wing_curvature","shape_complexity")
k_results <- lapply(traits_K, function(tr) {
  vals <- setNames(phylo_sp[[tr]], phylo_sp$tip)
  vals <- vals[!is.na(vals) & names(vals) %in% pruned$tip.label]
  if (length(unique(vals)) < 2) return(data.frame(trait=tr, n=length(vals), K=NA, p_value=NA))
  res <- tryCatch(phytools::phylosig(pruned, vals, method="K", test=TRUE, nsim=999), error=function(e) NULL)
  if (is.null(res)) return(data.frame(trait=tr, n=length(vals), K=NA, p_value=NA))
  data.frame(trait=tr, n=length(vals), K=round(res$K,6), p_value=round(res$P,6))
}) %>% bind_rows()
k_results <- bind_rows(k_results,
                       data.frame(trait="shapePC1", n=length(spc1), K=round(K1$K,6), p_value=round(K1$P,6)),
                       data.frame(trait="shapePC2", n=length(spc2), K=round(K2$K,6), p_value=round(K2$P,6))
)
write.csv(k_results, "clade_analyses_v4/results/phylogenetic_signal_K_v2.csv", row.names=FALSE)
cat("Blomberg K exported\n"); print(k_results)

p_all   <- make_phylo_v4("all")
p_nonp  <- make_phylo_v4("non_pterodact")
p_pt    <- make_phylo_v4("pterodact")
p_nonpf <- make_phylo_v4("non_pterodactyliformes")   # NEW — excludes Darwinoptera

ggsave("clade_analyses_v4/plots_PDF/CLADE_15A_phylo_non_pterodactyliformes.pdf",
       p_nonpf, width=11, height=8.5, dpi=300)
ggsave("clade_analyses_v4/plots_PDF/CLADE_15A_phylo_nonpterodact.pdf",
       p_nonp,  width=11, height=8.5, dpi=300)
ggsave("clade_analyses_v4/plots_PDF/CLADE_15B_phylo_pterodact.pdf",
       p_pt,    width=11, height=8.5, dpi=300)
ggsave("clade_analyses_v4/plots_PDF/CLADE_15_phylo_all.pdf",
       p_all,   width=11, height=8.5, dpi=300)
ggsave("clade_analyses_v4/plots_PDF/CLADE_15C_pareto_violin.pdf",
       p_violin, width=7, height=8, dpi=300)

# Combined: Non-Pterodactyliformes (a) + Pterodactyloidea (b) + violin (c)
combined_ab <- plot_grid(p_nonpf, p_pt, ncol=2,
                         labels=c("a","b"), label_size=11, label_fontface="bold")
ggsave("clade_analyses_v4/plots_PDF/CLADE_15_phylo_combined.pdf",
       combined_ab, width=22, height=9, dpi=300)
combined_abc <- plot_grid(combined_ab, p_violin, ncol=1, rel_heights=c(2.2,1),
                          labels=c("","c"), label_size=11, label_fontface="bold")
ggsave("clade_analyses_v4/plots_PDF/CLADE_15_phylo_abc.pdf",
       combined_abc, width=22, height=20, dpi=300)
cat("All phylomorphospace figures saved.\n")