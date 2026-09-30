# Measure disparity through time and build convex hulls (Part B, continues Msc_temp_3.1.R).
cat("\n================================================================================\n")
cat("TEMPORAL METRICS - PART B (disparity, EFA disparity, convex hulls)\n")
cat("================================================================================\n\n")

if (!exists("temporal_data")) stop("Source Msc_temp_3_1a.R first.")

# shape_scores lacks Midpoint/Period.Name - join them in from
# performance_data_clean by species before testing for them.
if (exists("shape_scores") && !all(c("Midpoint", "Period.Name") %in% colnames(shape_scores)) &&
    exists("performance_data_clean")) {
  shape_scores <- shape_scores %>%
    dplyr::left_join(
      performance_data_clean %>% dplyr::select(species, Midpoint, Period.Name) %>% dplyr::distinct(),
      by = "species")
  cat("Midpoint/Period.Name joined into shape_scores from performance_data_clean.\n")
}

######################################################################################
# 4. Disparity through time - PERFORMANCE PC plane (scale.=TRUE PCA, values
# naturally land around 1-10 after summing variances across a bootstrap
# sample - NOT comparable to the EFA panel in section 4b below).
######################################################################################
disparity_boot <- function(d, n_boot = 1000) {
  M <- cbind(d$PC1, d$PC2)
  M <- M[complete.cases(M), , drop = FALSE]
  if (nrow(M) < 4) return(c(sov = NA, sov_lo = NA, sov_hi = NA,
                            mpd = NA, mpd_lo = NA, mpd_hi = NA))
  sov <- function(m) sum(apply(m, 2, var))
  mpd <- function(m) mean(dist(m))

  bs <- replicate(n_boot, {
    idx <- sample(nrow(M), nrow(M), replace = TRUE)
    c(sov(M[idx, , drop = FALSE]), mpd(M[idx, , drop = FALSE]))
  })
  c(sov = sov(M), sov_lo = quantile(bs[1, ], 0.025), sov_hi = quantile(bs[1, ], 0.975),
    mpd = mpd(M), mpd_lo = quantile(bs[2, ], 0.025), mpd_hi = quantile(bs[2, ], 0.975))
}

disparity_metrics <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
  d <- temporal_data[temporal_data$Time_Bin == tb, ]
  if (nrow(d) == 0) return(NULL)
  cbind(data.frame(Time_Bin = tb, n = nrow(d)),
        as.data.frame(t(disparity_boot(d))))
}))
disparity_metrics$Time_Bin <- factor(disparity_metrics$Time_Bin,
                                     levels = levels(temporal_data$Time_Bin))
rownames(disparity_metrics) <- NULL

cat("Disparity by time bin (performance PC1-PC2 plane):\n")
print(disparity_metrics, digits = 4)
cat("\n")

######################################################################################
# 4b. Disparity through time - EFA SHAPE axes (shapePC1/shapePC2), NOT
# standardised with scale.=TRUE like the performance PCA above, so values
# stay naturally small. Requires shape_scores (species, Midpoint, Period.Name,
# shapePC1, shapePC2 - the join above adds Midpoint/Period.Name if missing).
######################################################################################
if (exists("shape_scores") &&
    all(c("shapePC1", "shapePC2", "Midpoint", "Period.Name") %in% colnames(shape_scores))) {

  efa_data <- shape_scores[!is.na(shape_scores$Midpoint) & !is.na(shape_scores$Period.Name), ]
  efa_data$Period_Combined <- efa_data$Period.Name
  efa_data$Period_Combined[efa_data$Period_Combined %in%
                             c("Early Jurassic", "Middle Jurassic")] <- "Early-Middle Jurassic"
  efa_data$Time_Bin <- factor(efa_data$Period_Combined,
                              levels = levels(temporal_data$Time_Bin))
  efa_data <- efa_data[!is.na(efa_data$Time_Bin), ]

  cat(sprintf("EFA shape-axis specimens with temporal data: %d\n", nrow(efa_data)))
  print(table(efa_data$Time_Bin))
  cat("\n")

  disparity_boot_efa <- function(d, n_boot = 1000) {
    M <- cbind(d$shapePC1, d$shapePC2)
    M <- M[complete.cases(M), , drop = FALSE]
    if (nrow(M) < 4) return(c(sov = NA, sov_lo = NA, sov_hi = NA,
                              mpd = NA, mpd_lo = NA, mpd_hi = NA))
    sov <- function(m) sum(apply(m, 2, var))
    mpd <- function(m) mean(dist(m))

    bs <- replicate(n_boot, {
      idx <- sample(nrow(M), nrow(M), replace = TRUE)
      c(sov(M[idx, , drop = FALSE]), mpd(M[idx, , drop = FALSE]))
    })
    c(sov = sov(M), sov_lo = quantile(bs[1, ], 0.025), sov_hi = quantile(bs[1, ], 0.975),
      mpd = mpd(M), mpd_lo = quantile(bs[2, ], 0.025), mpd_hi = quantile(bs[2, ], 0.975))
  }

  disparity_metrics_efa <- do.call(rbind, lapply(levels(efa_data$Time_Bin), function(tb) {
    d <- efa_data[efa_data$Time_Bin == tb, ]
    if (nrow(d) == 0) return(NULL)
    cbind(data.frame(Time_Bin = tb, n = nrow(d)),
          as.data.frame(t(disparity_boot_efa(d))))
  }))
  disparity_metrics_efa$Time_Bin <- factor(disparity_metrics_efa$Time_Bin,
                                           levels = levels(temporal_data$Time_Bin))
  rownames(disparity_metrics_efa) <- NULL

  cat("Disparity by time bin (EFA shapePC1-shapePC2 plane):\n")
  print(disparity_metrics_efa, digits = 4)
  cat("\n")
} else {
  disparity_metrics_efa <- NULL
  cat("shape_scores (with shapePC1/shapePC2/Midpoint/Period.Name) not found -\n")
  cat("EFA-axis disparity skipped.\n\n")
}

######################################################################################
# 5. Build the convex hull (outer boundary points) per clade and time bin,
# used to shade each clade's occupied region in the temporal figure.
######################################################################################
temporal_hulls_clade <- do.call(rbind, lapply(levels(temporal_data$Time_Bin), function(tb) {
  d <- temporal_data[temporal_data$Time_Bin == tb, ]
  do.call(rbind, lapply(unique(d$clade), function(cl) {
    dc <- d[d$clade == cl & is.finite(d$PC1) & is.finite(d$PC2), ]
    if (nrow(dc) < 3) return(NULL)
    h <- chull(dc$PC1, dc$PC2)
    data.frame(Time_Bin = tb, clade = cl, PC1 = dc$PC1[h], PC2 = dc$PC2[h])
  }))
}))
temporal_hulls_clade$Time_Bin <- factor(temporal_hulls_clade$Time_Bin, levels = levels(temporal_data$Time_Bin))
cat(sprintf("Clade hulls built: %d clade-bin combinations\n\n",nrow(unique(temporal_hulls_clade[, c("Time_Bin", "clade")]))))

cat("Available: disparity_metrics, disparity_metrics_efa, temporal_hulls_clade\n\n")
