# --- 1. Load the results table used for the preservation test -----------------
dat <- read.csv("output/results/performance_metrics_complete_CLADE.csv",
                stringsAsFactors = FALSE)

# --- 2. Turn depositional labels into preservation groups ---------------------
# Exact, mutually-exclusive matches against the categories present in the data.
# Order does not matter here because each pattern targets a distinct label,
# but we keep High -> Low to read top-down.
env_raw <- trimws(dat$Depositional.settings.paleoenvironment)

dat$Pres_Quality <- dplyr::case_when(
  grepl("Lagoon",                  env_raw, ignore.case = TRUE) ~ "High",
  grepl("Lacustrine:\\s*small",    env_raw, ignore.case = TRUE) ~ "High",
  grepl("Lacustrine:\\s*large",    env_raw, ignore.case = TRUE) ~ "Med-High",
  grepl("Fluviodeltaic",           env_raw, ignore.case = TRUE) ~ "Medium",
  grepl("Coastal|Shallow marine",  env_raw, ignore.case = TRUE) ~ "Medium",
  grepl("Alluvial",                env_raw, ignore.case = TRUE) ~ "Low",
  grepl("playa",                   env_raw, ignore.case = TRUE) ~ "Low",
  grepl("Aeolian",                 env_raw, ignore.case = TRUE) ~ "Low",
  TRUE ~ NA_character_   # anything unmatched is flagged, not silently "Low"
)

# Sanity check: report any unmatched depositional settings before dropping them
unmatched <- unique(env_raw[is.na(dat$Pres_Quality) & env_raw != ""])
if (length(unmatched) > 0) {
  warning("Unmatched depositional settings (set to NA): ",
          paste(unmatched, collapse = " | "))
}

# Keep the groups in a fixed order before plotting and testing.
dat$Pres_Quality <- factor(dat$Pres_Quality,
                           levels = c("Low", "Medium", "Med-High", "High"))

dat <- dat %>% filter(!is.na(Pres_Quality) & !is.na(PC2))

# Quick cross-tab so you can confirm the mapping / sample sizes
cat("\nSpecimens per preservation-quality category:\n")
print(table(dat$Pres_Quality))
cat("\nDepositional setting x preservation quality:\n")
print(table(trimws(dat$Depositional.settings.paleoenvironment), dat$Pres_Quality))

# --- 3. Group means (for annotations) ----------------------------------------
group_means <- dat %>%
  group_by(Pres_Quality) %>%
  summarise(mean_PC2 = mean(PC2, na.rm = TRUE),
            n = n(), .groups = "drop")

# --- 4. Kruskal-Wallis -------------------------------------------------------
kw <- kruskal.test(PC2 ~ Pres_Quality, data = dat)
print(kw)

# --- 5. Plot PC2 for each preservation group ----------------------------------
y_lab_pos <- max(dat$PC2, na.rm = TRUE) * 1.05

p1 <- ggplot(dat, aes(x = Pres_Quality, y = PC2, fill = Pres_Quality)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.18, alpha = 0.6, size = 2, aes(color = Pres_Quality)) +
  scale_fill_viridis(discrete = TRUE, end = 0.8) +
  scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
  ylab("PC2 score (wing complexity/curvature)") +
  xlab("Preservation quality (depositional environment)") +
  ggtitle("Taphonomic sensitivity: PC2 by preservation quality") +
  geom_text(data = group_means,
            aes(x = Pres_Quality, y = y_lab_pos,
                label = paste0("\u03bc=", round(mean_PC2, 2), "\n(n=", n, ")")),
            size = 3.5, inherit.aes = FALSE) +
  theme_minimal(base_size = 14) +
  guides(color = "none", fill = "none")

# --- 6. Test the simple rank relationship -------------------------------------
# Code Low..High as 1..4 to match the factor level order above.
dat$Pres_Score <- as.numeric(dat$Pres_Quality)

corsimple <- cor.test(dat$PC2, dat$Pres_Score, method = "spearman")
cat(sprintf("\nRaw correlation: rho=%.3f, p=%.3g\n",
            corsimple$estimate, corsimple$p.value))

p2 <- ggplot(dat, aes(x = Pres_Quality, y = PC2)) +
  geom_jitter(width = 0.18, size = 2, alpha = 0.7, aes(color = Pres_Quality)) +
  stat_smooth(method = "lm", aes(group = 1), size = 1.1, se = TRUE, color = "black") +
  scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
  ylab("PC2 score") +
  xlab("Preservation quality") +
  annotate("text", x = 1, y = max(dat$PC2, na.rm = TRUE),
           label = sprintf("rho=%.2f, p=%.2g", corsimple$estimate, corsimple$p.value),
           size = 5, hjust = 0) +
  theme_minimal(base_size = 13) + guides(color = "none")

# --- 7. Linear model: control for clade (partial correlation) ----------------
mod <- lm(PC2 ~ Pres_Score + clade, data = dat)
cat("\n--- lm(PC2 ~ Pres_Score + clade) ---\n")
print(summary(mod))

dat$PC2_clade_resid      <- residuals(lm(PC2 ~ clade, data = dat))
dat$PresScore_clade_resid <- residuals(lm(Pres_Score ~ clade, data = dat))

corresid <- cor.test(dat$PC2_clade_resid, dat$PresScore_clade_resid,
                     method = "spearman")
cat(sprintf("\nPartial correlation after clade: rho=%.3f, p=%.3g\n",
            corresid$estimate, corresid$p.value))

p3 <- ggplot(dat, aes(x = PresScore_clade_resid, y = PC2_clade_resid,
                      color = Pres_Quality)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_viridis(discrete = TRUE, end = 0.8, option = "D") +
  xlab("Preservation quality (residuals after clade)") +
  ylab("PC2 (residuals after clade)") +
  annotate("text",
           x = min(dat$PresScore_clade_resid, na.rm = TRUE),
           y = max(dat$PC2_clade_resid, na.rm = TRUE),
           label = sprintf("rho=%.2f, p=%.2g", corresid$estimate, corresid$p.value),
           hjust = 0, size = 5) +
  theme_minimal(base_size = 13) + guides(color = "none")

# --- 8. Heatmap: mean PC2 by clade & preservation ----------------------------
hm <- dat %>%
  group_by(clade, Pres_Quality) %>%
  summarise(mean_PC2 = mean(PC2, na.rm = TRUE), n = n(), .groups = "drop")

hm$label <- with(hm, sprintf("%.1f\n(n=%d)", mean_PC2, n))

# keep all 4 preservation levels on the x-axis even where a clade has no data
hm$Pres_Quality <- factor(hm$Pres_Quality,
                          levels = c("Low", "Medium", "Med-High", "High"))
hm$clade <- factor(hm$clade, levels = rev(sort(unique(dat$clade))))

p4 <- ggplot(hm, aes(x = Pres_Quality, y = clade, fill = mean_PC2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "#1a9850", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Mean PC2") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_size = 12) +
  labs(x = "Preservation quality", y = "Clade") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

# --- 9. Output ---------------------------------------------------------------
fig_box <- ggarrange(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"))
ggsave("output/plots_PDF/taphonomic_sensitivity_pc2_box_corrs.pdf", fig_box,
       width = 14, height = 5, dpi = 300)

ggsave("output/plots_PDF/taphonomic_sensitivity_heatmap.pdf", p4, width = 8, height = 5, dpi = 300)

# --- Optional: summaries for report ------------------------------------------
# Show a fixed number of decimals so every row lines up.
hm_display <- hm
hm_display$mean_PC2 <- format_fixed(hm$mean_PC2, 3)
write.csv(hm_display, "output/results/supplementals/taphonomic_PC2_by_clade_prescat.csv", row.names = FALSE)
