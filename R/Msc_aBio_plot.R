cat("\n================================================================================\n")
cat("MORPHOSPACE OCCUPATION BY GROUP\n")
cat("================================================================================\n\n")

# --------------------------------------------------------------------------------
# 1. Diet combination — the category the third figure needs
# --------------------------------------------------------------------------------
if (all(c("Diet.1", "Diet.2") %in% colnames(performance_data_clean))) {
  d1 <- trimws(as.character(performance_data_clean$Diet.1))
  d2 <- trimws(as.character(performance_data_clean$Diet.2))
  performance_data_clean$Diet_combined <- ifelse(
    is.na(d2) | d2 == "" | d2 == d1, d1, paste(d1, d2, sep = " + "))
  cat("Diet combinations found:", length(unique(na.omit(
    performance_data_clean$Diet_combined))), "\n\n")
}

GROUPINGS <- list(
  list(col = "clade", palette = "Set3", fixed = CLADE_COLS,
       title = "Morphospace occupation by clade",
       file  = "CLADE_10_morphospace_clade_hulls"),
  list(col = "Environment", palette = "Set2", fixed = ENV_COLS,
       title = "Morphospace occupation by environment",
       file  = "PALEV_10_morphospace_environment_hulls"),
  list(col = "Diet.1", palette = "Dark2", fixed = DIET_COLS,
       title = "Morphospace occupation by primary diet",
       file  = "DIET_10_morphospace_diet_hulls"),
  list(col = "Diet_combined", palette = "Paired", fixed = c(DIETCOMBO_COLS, DIET_COLS),
       title = "Morphospace occupation by diet combination",
       file  = "DIET_11_morphospace_diet_combination_hulls"),
  list(col = "Depositional.settings.paleoenvironment", palette = "Set2", fixed = DEPOSITIONAL_COLS,
       title = "Morphospace Occupation by Depositional Setting",
       file  = "PALEV_12_morphospace_depositional_convex_hulls")
)

# --------------------------------------------------------------------------------
# 2. One builder for every grouping
# --------------------------------------------------------------------------------
make_group_plot <- function(group_col, palette_name, title_text, fixed_cols = NULL) {
  
  if (!group_col %in% colnames(performance_data_clean)) {
    cat("SKIPPED", group_col, "- column not found\n"); return(NULL)
  }
  
  d <- performance_data_clean[!is.na(performance_data_clean[[group_col]]) &
                                is.finite(performance_data_clean$PC1) &
                                is.finite(performance_data_clean$PC2), ]
  d$.grp <- as.character(d[[group_col]])
  d <- d[d$.grp != "" & !is.na(d$.grp), ]
  if (nrow(d) < 5) { cat("SKIPPED", group_col, "- too few specimens\n"); return(NULL) }
  
  groups <- sort(unique(d$.grp)); n_g <- length(groups)
  
  # Fixed colour codes take priority over the auto-generated Brewer palette,
  # matched case/whitespace-insensitively so labelling differences in the
  # data don't silently fall through to the fallback palette.
  if (!is.null(fixed_cols)) {
    cols <- match_fixed_colour(groups, fixed_cols)
    missing_cl <- groups[is.na(cols)]
    if (length(missing_cl) > 0) {
      cat("  NOTE: no fixed colour entry for:", paste(missing_cl, collapse = ", "),
          "- falling back to auto palette for these.\n")
      fallback <- if (n_g <= 8) {
        brewer.pal(max(3, n_g), palette_name)[seq_len(max(3, n_g))][seq_len(n_g)]
      } else {
        colorRampPalette(brewer.pal(min(12, max(3, n_g)), palette_name))(n_g)
      }
      names(fallback) <- groups
      cols[missing_cl] <- fallback[missing_cl]
    }
    # DIAGNOSTIC: print the exact colour used for each group.
    cat("  Colour mapping used:\n"); print(cols)
  } else {
    cols <- if (n_g <= 8) {
      brewer.pal(max(3, n_g), palette_name)[seq_len(max(3, n_g))][seq_len(n_g)]
    } else {
      colorRampPalette(brewer.pal(min(12, max(3, n_g)), palette_name))(n_g)
    }
    names(cols) <- groups
  }
  
  # Hulls for n >= 3, segments for n == 2, nothing for n == 1
  hulls <- do.call(rbind, lapply(groups, function(g) {
    dg <- d[d$.grp == g, ]
    if (nrow(dg) < 3) return(NULL)
    h <- chull(dg$PC1, dg$PC2); h <- c(h, h[1])
    data.frame(PC1 = dg$PC1[h], PC2 = dg$PC2[h], .grp = g)
  }))
  
  pairs2 <- do.call(rbind, lapply(groups, function(g) {
    dg <- d[d$.grp == g, ]
    if (nrow(dg) != 2) return(NULL)
    data.frame(PC1 = dg$PC1, PC2 = dg$PC2, .grp = g)
  }))
  
  n_h <- if (is.null(hulls))  0 else length(unique(hulls$.grp))
  n_p <- if (is.null(pairs2)) 0 else length(unique(pairs2$.grp))
  cat(sprintf("  %s: %d groups | %d hulls | %d segments | %d single points\n",
              group_col, n_g, n_h, n_p, n_g - n_h - n_p))
  
  p <- ggplot() +
    geom_tile(data = stress_grid, aes(x = PC1, y = PC2, fill = stress),
              alpha = 0.5) +
    scale_fill_stress(palette = "dark", name = surf_label)
  
  if (!is.null(wing_bg)) {
    p <- p + geom_polygon(data = wing_bg, aes(x = PC1, y = PC2, group = grid_id),
                          fill = "#2a2a2a", colour = "#1a1a1a",
                          alpha = 0.16, linewidth = 0.25)
  }
  
  p <- p + ggnewscale::new_scale_fill()
  
  if (!is.null(hulls)) {
    p <- p +
      geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = .grp, fill = .grp),
                   alpha = 0.26, colour = NA) +
      geom_path(data = hulls, aes(x = PC1, y = PC2, group = .grp, colour = .grp),
                linewidth = 1.2, alpha = 0.9)
  }
  if (!is.null(pairs2)) {
    p <- p + geom_line(data = pairs2, aes(x = PC1, y = PC2, group = .grp, colour = .grp), linewidth = 1.0, alpha = 0.8, linetype = "dashed")
  }
  
  p + geom_point(data = d, aes(x = PC1, y = PC2, colour = .grp), size = 2.8, alpha = 0.85) + scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(title = title_text, subtitle = paste0("Surface: ", gsub("\n", " ", surf_label),". Dashed lines: groups of two specimens."),
         x = pc_lab(1), y = pc_lab(2)) + theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0, face = "bold", size = 16),
                                                                               plot.subtitle = element_text(hjust = 0, size = 9, colour = "grey35"),legend.position = "right",
                                                                               legend.text = element_text(size = 8),panel.grid.major = element_line(colour = "grey90"),
                                                                               panel.background = element_rect(fill = "white", colour = NA),plot.background  = element_rect(fill = "white", colour = NA),
                                                                               panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
}

# --------------------------------------------------------------------------------
# 3. Build and write
# --------------------------------------------------------------------------------
for (g in GROUPINGS) {
  p <- make_group_plot(g$col, g$palette, g$title, g$fixed)
  if (is.null(p)) next
  ggsave(sprintf("clade_analyses_v4/plots/biomechanic_plot/%s.png", g$file),p, width = 16, height = 11, dpi = 300)
  ggsave(sprintf("clade_analyses_v4/plots_PDF/biomechanic_plot/%s.pdf", g$file), p, width = 16, height = 11)
  cat("Written:", g$file, "\n")
}
cat("\nGroup morphospace figures complete.\n\n")