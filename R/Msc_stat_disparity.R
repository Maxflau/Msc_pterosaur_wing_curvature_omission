
if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")
if (!all(c("PC1", "PC2") %in% colnames(performance_data_clean))) {
  stop("PC1/PC2 not found - source Msc_impossible_region_v3.R first.")
}

#' Disparity of one grouping in the PC1-PC2 plane
disparity_by <- function(d, group_col) {
  
  d <- d[!is.na(d[[group_col]]) & d[[group_col]] != "" &
           is.finite(d$PC1) & is.finite(d$PC2), ]
  if (nrow(d) < MIN_N) return(NULL)
  
  out <- do.call(rbind, lapply(unique(d[[group_col]]), function(lv) {
    s <- d[d[[group_col]] == lv, ]
    if (nrow(s) < MIN_N) return(NULL)
    data.frame(group = lv, n = nrow(s),
               mean_PC1 = round(mean(s$PC1), 4),
               mean_PC2 = round(mean(s$PC2), 4),
               var_PC1  = round(var(s$PC1), 4),
               var_PC2  = round(var(s$PC2), 4),
               sum_of_variance = round(var(s$PC1) + var(s$PC2), 4),
               mean_pairwise_distance = round(mean(dist(cbind(s$PC1, s$PC2))), 4),
               stringsAsFactors = FALSE)
  }))
  
  if (is.null(out)) return(NULL)
  out <- out[order(-out$sum_of_variance), ]
  names(out)[1] <- group_col
  out
}

# --------------------------------------------------------------------------------
# 1. Taxonomic and ecological groupings
# --------------------------------------------------------------------------------
DISP_GROUPS <- intersect(c("clade", "Order", "Diet.1", "Diet_combined",
                           "Depositional.settings.paleoenvironment", "Environment"),
                         colnames(performance_data_clean))

for (g in DISP_GROUPS) {
  disp <- disparity_by(performance_data_clean, g)
  cat(sprintf("Disparity by %s: %d groups with n >= %d\n",
              g, if (is.null(disp)) 0 else nrow(disp), MIN_N))
  write_supp(disp, paste0("S4_disparity_by_",
                          tolower(gsub("\\.", "_", g))))
}

# Diet_combined is built in Msc_stats_00_setup.R and is already in DISP_GROUPS
# above, so it needs no separate block here.

cat("\nDisparity is computed on the performance PCs, so it reflects functional\n")
cat("spread rather than outline shape. For shape disparity use the EFA harmonics.\n\n")