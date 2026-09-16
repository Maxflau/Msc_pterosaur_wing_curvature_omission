# --- summary_pterosaur_env_performance_clean.R ---

# If using after the complete_pterosaur_phylomorpho2.R workflow
# (Assumes: performance_data_clean exists!)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(grid)
library(viridis)
library(ggnewscale)
library(patchwork)

# ------ 1. Source data: performance_data_clean, patch in deposition/env if missing

# It should contain 'Environment' and 'Depositional.settings.paleoenvironment'
# If you're running this after the big pipeline script, it's likely these are present.
# If not, join them in from your CSV:

# ------ 1. Patch in deposition/env from CSV if missing -----------------------
if(!("Environment" %in% colnames(performance_data_clean)) ||
   !("Depositional.settings.paleoenvironment" %in% colnames(performance_data_clean))) {
  master_df <- read.csv("data/pteros_main_data.csv", sep=";", stringsAsFactors=FALSE, check.names=FALSE)
  colnames(master_df) <- trimws(colnames(master_df))
  patch <- master_df %>% transmute(
    species = trimws(as.character(`SPECIES (Pegas)`)),
    Period.Name = trimws(as.character(`Period Name`)),
    Environment = trimws(as.character(Environment)),
    Depositional.settings.paleoenvironment = trimws(as.character(`Depositional settings/paleoenvironment`)))
  performance_data_clean <- performance_data_clean %>% left_join(patch, by = "species")
}

performance_data_clean$Period.Name <- trimws(as.character(performance_data_clean$Period.Name))
performance_data_clean$Environment <- trimws(as.character(performance_data_clean$Environment))
performance_data_clean$Depositional.settings.paleoenvironment <-
  trimws(as.character(performance_data_clean$Depositional.settings.paleoenvironment))

# Collapse Early + Middle Jurassic
performance_data_clean$Period_collapsed <- dplyr::recode(performance_data_clean$Period.Name,
                                                         "Early Jurassic" = "Early+Middle Jurassic", "Middle Jurassic" = "Early+Middle Jurassic",
                                                         .default = performance_data_clean$Period.Name)
period_levels <- c("Late Triassic", "Early+Middle Jurassic", "Late Jurassic",
                   "Early Cretaceous", "Late Cretaceous")
performance_data_clean$Period_collapsed <- factor(performance_data_clean$Period_collapsed, levels = period_levels)

# ------ 2. Tabulate by period ------------------------------------------------
tab_by_period <- function(df, col) {
  df %>% filter(!is.na(.data[[col]]), .data[[col]] != "", !is.na(Period_collapsed)) %>%
    group_by(Period_collapsed, cat = .data[[col]]) %>% summarise(N = n(), .groups="drop") %>%
    group_by(Period_collapsed) %>% mutate(prop = N / sum(N), total = sum(N)) %>% ungroup()
}
tab_depo <- tab_by_period(performance_data_clean, "Depositional.settings.paleoenvironment")
tab_env  <- tab_by_period(performance_data_clean, "Environment")

# ------ 3. Palettes (one colour per category present in the data) ------------
# 8 depositional settings: Aeolian, Alluvial plain, Coastal/Shallow marine,
# Fluviodeltaic, Lacustrine large/small, Lagoonal deposit, playa
depo_colors <- c(
  "Fluviodeltaic"="#117a65", "Alluvial plain"="#f39c12", "Coastal/Shallow marine"="#48c9b0",
  "Lagoonal deposit"="#aed6f1", "Lacustrine: small lake"="#5499c7", "Lacustrine: large lake"="#154360",
  "Aeolian"="#d35400", "playa"="#fcf3cf")
# 8 habitats: Archipelago, Coastal env, Desert, Floodplain, Fluvial env,
# Forested wetland env, Marine env, Semi-arid floodplain
env_colors <- c(
  "Coastal environment"="#23b0b5", "Forested wetland environment"="#078100",
  "Semi-arid floodplain"="#ac612a", "Fluvial environment"="#e67e22", "Floodplain"="#73c6b6",
  "Archipelago"="#bb8fce", "Marine environment"="#00bcd4", "Desert"="#e1b12c")

# Safety check: warn if any category in the data has no colour (would drop from legend)
check_palette <- function(levels_in_data, palette, label) {
  missing <- setdiff(unique(levels_in_data), names(palette))
  if (length(missing) > 0) warning(label, " - no colour for: ", paste(missing, collapse=" | "))
  else cat(sprintf("  %s: all %d categories coloured.\n", label, length(unique(levels_in_data))))
}
check_palette(tab_depo$cat, depo_colors, "Depositional setting")
check_palette(tab_env$cat,  env_colors,  "Habitat")

# ------ 4. Barplots ----------------------------------------------------------
make_bar <- function(tab, palette, fill_name, title, lab_thresh) {
  ggplot(tab, aes(x=Period_collapsed, y=prop, fill=cat)) +
    geom_bar(stat="identity", width=0.7, color="black") +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    scale_fill_manual(values=palette, name=fill_name, drop=FALSE) +
    geom_text(aes(label=ifelse(prop>lab_thresh, scales::percent(prop, accuracy=1), "")),
              position=position_stack(vjust=0.5), size=3, color="black") +
    geom_text(data=tab %>% group_by(Period_collapsed) %>% filter(row_number()==1),
              aes(y=1, label=paste0("n=", total)), vjust=-0.5, hjust=1, size=4, color="black") +
    labs(x="Period", y="Specimens (%)", title=title) +
    theme_minimal(base_size=13) + theme(axis.text.x=element_text(angle=30, hjust=1))
}
p_depo <- make_bar(tab_depo, depo_colors, "Depositional setting", "Depositional setting by period", 0.06)
p_env  <- make_bar(tab_env,  env_colors,  "Habitat", "Habitat (macroenvironment) by period", 0.07)

# ------ 5. Relationship: habitat vs depositional setting ---------------------
rel <- performance_data_clean %>%
  filter(!is.na(Environment), Environment != "",
         !is.na(Depositional.settings.paleoenvironment),
         Depositional.settings.paleoenvironment != "")
ct <- table(rel$Environment, rel$Depositional.settings.paleoenvironment)
cat("\nHabitat x depositional setting contingency table:\n"); print(ct)
# Chi-square (with Monte-Carlo p, robust to sparse cells) + Cramer's V effect size
chi <- chisq.test(ct, simulate.p.value = TRUE, B = 10000)
cramers_v <- sqrt(as.numeric(chi$statistic) / (sum(ct) * (min(dim(ct)) - 1)))
cat(sprintf("\nChi-square (MC): X2=%.1f, p=%.4g | Cramer's V=%.3f\n",
            chi$statistic, chi$p.value, cramers_v))

# ------ 6. Combined plot + save ----------------------------------------------
(p_depo | p_env) + plot_annotation(title = "Summary of pterosaur environments through time")
ggsave("output/plots/Taphono_plot.png", p_depo + p_env, width=14, height=8, dpi=300)
ggsave("output/plots_PDF/Taphono_plot.pdf", p_depo + p_env, width=14, height=8)