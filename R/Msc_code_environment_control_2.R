# Compare environments, depositional settings, and performance metrics.

# ------ 1. Source data: performance_data_clean, patch in deposition/env if missing

# ------ 1. Patch in deposition/env from CSV if missing -----------------------
# If needed, copy missing environment labels from the master CSV.
if(!all(c("Environment","Depositional.settings.paleoenvironment") %in% colnames(performance_data_clean))) {
  master_df <- read.csv("data/pteros_main_data.csv", sep=";", stringsAsFactors=FALSE, check.names=FALSE)
  colnames(master_df) <- trimws(colnames(master_df))
  patch <- master_df %>% transmute(species=trimws(`SPECIES (Pegas)`), Period.Name=trimws(`Period Name`),
                                   Environment=trimws(Environment), Depositional.settings.paleoenvironment=trimws(`Depositional settings/paleoenvironment`))
  performance_data_clean <- left_join(performance_data_clean, patch, by="species")
}
for(c in c("Period.Name","Environment","Depositional.settings.paleoenvironment"))
  performance_data_clean[[c]] <- trimws(as.character(performance_data_clean[[c]]))

# Merge the two Jurassic bins so the time categories stay consistent.
performance_data_clean$Period_collapsed <- factor(dplyr::recode(performance_data_clean$Period.Name,
                                                                "Early Jurassic"="Early+Middle Jurassic", "Middle Jurassic"="Early+Middle Jurassic",
                                                                .default=performance_data_clean$Period.Name),
                                                  levels=c("Late Triassic","Early+Middle Jurassic","Late Jurassic","Early Cretaceous","Late Cretaceous"))

# ------ 2. Tabulate by period ------------------------------------------------
tab_by_period <- function(df, col) df %>%
  filter(!is.na(.data[[col]]), .data[[col]]!="", !is.na(Period_collapsed)) %>%
  group_by(Period_collapsed, cat=.data[[col]]) %>% summarise(N=n(), .groups="drop") %>%
  group_by(Period_collapsed) %>% mutate(prop=N/sum(N), total=sum(N)) %>% ungroup()
tab_depo <- tab_by_period(performance_data_clean, "Depositional.settings.paleoenvironment")
tab_env  <- tab_by_period(performance_data_clean, "Environment")

# ------ 3. Set the colours used in the figures --------------------------------
depo_colors <- c("Fluviodeltaic"="#117a65",
                 "Alluvial plain"="#f39c12",
                 "Coastal/Shallow marine"="#48c9b0",
                 "Lagoonal deposit"="#aed6f1",
                 "Lacustrine: small lake"="#5499c7",
                 "Lacustrine: large lake"="#154360",
                 "Aeolian"="#d35400",
                 "playa"="#fcf3cf")
env_colors <- c("Coastal environment"="#23b0b5","Forested wetland environment"="#078100","Semi-arid floodplain"="#ac612a",
                "Fluvial environment"="#e67e22","Floodplain"="#73c6b6","Archipelago"="#bb8fce","Marine environment"="#00bcd4","Desert"="#e1b12c")
check_palette <- function(lv, pal, lab) { miss <- setdiff(unique(lv), names(pal))
if(length(miss)) warning(lab," - no colour for: ", paste(miss, collapse=" | "))
else cat(sprintf("  %s: all %d categories coloured.\n", lab, length(unique(lv)))) }
check_palette(tab_depo$cat, depo_colors, "Depositional setting"); check_palette(tab_env$cat, env_colors, "Habitat")

# ------ 4. Barplots ----------------------------------------------------------
make_bar <- function(tab, palette, fill_name, title, thr) ggplot(tab, aes(Period_collapsed, prop, fill=cat)) +
  geom_bar(stat="identity", width=0.7, color="black") +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  scale_fill_manual(values=palette, name=fill_name, drop=FALSE) +
  geom_text(aes(label=ifelse(prop>thr, scales::percent(prop, accuracy=1), "")), position=position_stack(vjust=0.5), size=3) +
  geom_text(data=tab %>% group_by(Period_collapsed) %>% filter(row_number()==1),
            aes(y=1, label=paste0("n=", total)), vjust=-0.5, hjust=1, size=4) +
  labs(x="Period", y="Specimens (%)", title=title) +
  theme_minimal(base_size=13) + theme(axis.text.x=element_text(angle=30, hjust=1))
p_depo <- make_bar(tab_depo, depo_colors, "Depositional setting", "Depositional setting by period", 0.06)
p_env  <- make_bar(tab_env,  env_colors,  "Habitat", "Habitat (macroenvironment) by period", 0.07)
(p_depo | p_env) + plot_annotation(title="Summary of pterosaur environments through time")
ggsave("output/plots/Taphono_plot.png", p_depo + p_env, width=14, height=8, dpi=300)
ggsave("output/plots_PDF/Taphono_plot.pdf", p_depo + p_env, width=14, height=8)

# ------ 5. Check whether habitat and depositional setting overlap -------------
rel <- performance_data_clean %>% filter(!is.na(Environment), Environment!="",
                                         !is.na(Depositional.settings.paleoenvironment), Depositional.settings.paleoenvironment!="")
ct <- table(rel$Environment, rel$Depositional.settings.paleoenvironment)
chi <- chisq.test(ct, simulate.p.value=TRUE, B=10000)
cramers_v <- sqrt(as.numeric(chi$statistic) / (sum(ct) * (min(dim(ct)) - 1)))
cat(sprintf("\nHabitat x deposition: X2=%.1f, p=%.4g | Cramer's V=%.3f\n", chi$statistic, chi$p.value, cramers_v))

# ------ 6. Test biomechanics against setting and time -------------------------
# Use the current column names so missing legacy names do not stop the script.
# second_moment -> r2_hat (was never a real column, so dd[["second_moment"]]
# was NULL and lm()/model.frame() could not build the model).
# wing_loading -> wing_loading_ratio (the current size-corrected convention).
# mean_optimality is kept only as a fallback if the legacy score is what this
# session happens to contain instead of pareto_rank_ratio.
metrics_wanted <- c("aspect_ratio", "wing_loading_ratio", "von_mises_stress",
                    "r2_hat", "pareto_rank_ratio", "mean_optimality")
metrics <- intersect(metrics_wanted, colnames(performance_data_clean))
missing_metrics <- setdiff(metrics_wanted, metrics)
if (length(missing_metrics) > 0) {
  cat("NOTE: skipping metrics not found in performance_data_clean:",
      paste(missing_metrics, collapse = ", "), "\n")
}
if (length(metrics) == 0) {
  stop("None of the wanted metrics are present - check performance_data_clean.")
}
cat("Metrics used:", paste(metrics, collapse = ", "), "\n\n")

dd <- rel
# Test group differences, time trends, and a combined linear model.
results <- lapply(metrics, function(m) {
  lm_fit <- lm(dd[[m]] ~ dd$Midpoint + factor(dd$Environment) + factor(dd$Depositional.settings.paleoenvironment))
  data.frame(metric=m,
             KW_habitat_p=kruskal.test(dd[[m]] ~ factor(dd$Environment))$p.value,
             KW_depo_p=kruskal.test(dd[[m]] ~ factor(dd$Depositional.settings.paleoenvironment))$p.value,
             KW_period_p=kruskal.test(dd[[m]] ~ dd$Period_collapsed)$p.value,
             time_rho=suppressWarnings(cor.test(dd[[m]], dd$Midpoint, method="spearman")$estimate),
             time_p=suppressWarnings(cor.test(dd[[m]], dd$Midpoint, method="spearman")$p.value),
             lm_adjR2=summary(lm_fit)$adj.r.squared,
             lm_time_p=summary(lm_fit)$coefficients["dd$Midpoint","Pr(>|t|)"])
}) %>% bind_rows() %>% mutate(across(where(is.numeric), ~round(.x, 4)))
# Show a fixed number of decimals so every row lines up.
results_display <- results %>% mutate(across(where(is.numeric), ~format_fixed(.x, 4)))
cat("\n=== Biomechanics vs environment / deposition / time ===\n"); print(results_display)
write.csv(results_display, "output/results/biomech_env_time_tests.csv", row.names=FALSE)

# Plot each metric through time, coloured by habitat.
trend_long <- dd %>% select(Midpoint, Environment, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to="metric", values_to="value")
p_trend <- ggplot(trend_long, aes(Midpoint, value)) +
  geom_point(aes(color=Environment), size=2, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.8) +
  scale_color_manual(values=env_colors, name="Habitat", drop=FALSE) +
  scale_x_reverse() + facet_wrap(~metric, scales="free_y") +
  labs(x="Age (Ma)", y="Metric value", title="Biomechanical performance through time by habitat") +
  theme_minimal(base_size=12)
ggsave("output/plots/Biomech_env_time.png", p_trend, width=12, height=8, dpi=300)
ggsave("output/plots_PDF/Biomech_env_time.pdf", p_trend, width=12, height=8)
