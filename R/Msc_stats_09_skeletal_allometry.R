cat("\n--- S12/S13: SKELETAL ALLOMETRY ---\n\n")

# Test how each bone changes as overall skeletal size increases.
if (!exists("POSTCRANIAL")) stop("Source Msc_stats_06a_postcranial.R first.")

MIN_FIT <- 10   # fewer points than this gives an unstable slope

#' Log-log regression of one element on skeletal size, tested against isometry.
#' Both variables are lengths, so the isometric slope is 1.
fit_element <- function(el, d, size_var = "SkeletalSize", expected = 1) {
  x <- d[[size_var]]; y <- d[[el]]
  ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  if (sum(ok) < MIN_FIT) return(NULL)

  m  <- lm(log(y[ok]) ~ log(x[ok]))
  cf <- summary(m)$coefficients
  slope <- cf[2, 1]; se <- cf[2, 2]

  # Confidence interval on the slope, so isometry can be judged by eye too
  ci <- slope + c(-1, 1) * qt(0.975, sum(ok) - 2) * se
  t_iso <- (slope - expected) / se
  p_iso <- 2 * pt(abs(t_iso), df = sum(ok) - 2, lower.tail = FALSE)

  data.frame(element = el, n = sum(ok),
             slope = round(slope, 4), se = round(se, 4),
             ci_lower = round(ci[1], 4), ci_upper = round(ci[2], 4),
             r_squared = round(summary(m)$r.squared, 4),
             expected_slope = expected,
             p_vs_isometry = signif(p_iso, 4),
             stringsAsFactors = FALSE)
}

# --------------------------------------------------------------------------------
# 1. Global allometry of each element
# --------------------------------------------------------------------------------
# Fit one size relationship for each measured bone across the full sample.
allo <- do.call(rbind, lapply(GM_ELEMENTS, fit_element, d = POSTCRANIAL))
if (is.null(allo)) stop("No element fit succeeded - check SkeletalSize.")

allo$p_adjusted_BH <- signif(p.adjust(allo$p_vs_isometry, method = "BH"), 4)
allo$allometry <- ifelse(allo$p_adjusted_BH >= 0.05, "isometric",
                         ifelse(allo$slope > allo$expected_slope,
                                "positive allometry", "negative allometry"))
allo <- allo[order(-allo$slope), ]

cat("Element allometry against the geometric-mean skeletal size:\n")
# Show a fixed number of decimals so every row lines up.
allo_display <- allo
allo_display$slope         <- format_fixed(allo$slope, 4)
allo_display$ci_lower      <- format_fixed(allo$ci_lower, 4)
allo_display$ci_upper      <- format_fixed(allo$ci_upper, 4)
allo_display$r_squared     <- format_fixed(allo$r_squared, 4)
allo_display$p_adjusted_BH <- format_fixed(allo$p_adjusted_BH, 4)
print(allo_display[, c("element", "n", "slope", "ci_lower", "ci_upper", "r_squared",
               "p_adjusted_BH", "allometry")], row.names = FALSE)
write_supp(allo_display, "S12_skeletal_allometry")

cat("\nSlope 1 = the element keeps its proportion as the animal grows.\n")
cat("Above 1 = it grows faster than the skeleton as a whole; below 1, slower.\n")
cat("This is where a disproportionately elongated wing finger would show up.\n\n")

# --------------------------------------------------------------------------------
# 2. Forelimb against hindlimb
# --------------------------------------------------------------------------------
# The biologically loaded contrast in pterosaurs: does the flight apparatus
# scale differently from the walking apparatus?
FORE <- intersect(c("Humerus", "Ulna", "McIV", "WingPh1", "WingPh2",
                    "WingPh3", "WingPh4"), GM_ELEMENTS)
HIND <- intersect(c("Femur", "Tibia", "MtIII"), GM_ELEMENTS)

if (length(FORE) >= 3 && length(HIND) >= 2) {
  POSTCRANIAL$ForelimbGM <- exp(rowMeans(log(POSTCRANIAL[, FORE]), na.rm = TRUE))
  POSTCRANIAL$HindlimbGM <- exp(rowMeans(log(POSTCRANIAL[, HIND]), na.rm = TRUE))

  ok <- is.finite(POSTCRANIAL$ForelimbGM) & is.finite(POSTCRANIAL$HindlimbGM) &
    POSTCRANIAL$ForelimbGM > 0 & POSTCRANIAL$HindlimbGM > 0

  m <- lm(log(POSTCRANIAL$ForelimbGM[ok]) ~ log(POSTCRANIAL$HindlimbGM[ok]))
  cf <- summary(m)$coefficients
  p_iso <- 2 * pt(abs((cf[2, 1] - 1) / cf[2, 2]), df = sum(ok) - 2, lower.tail = FALSE)

  cat(sprintf("Forelimb on hindlimb: slope = %.3f (SE %.3f), R² = %.3f, n = %d\n",
              cf[2, 1], cf[2, 2], summary(m)$r.squared, sum(ok)))
  cat(sprintf("Against isometry (slope 1): p = %.4g -> %s\n\n", p_iso,
              ifelse(p_iso >= 0.05, "isometric",
                     ifelse(cf[2, 1] > 1, "forelimb grows faster",
                            "hindlimb grows faster"))))

  POSTCRANIAL$ForeHindRatio <- POSTCRANIAL$ForelimbGM / POSTCRANIAL$HindlimbGM
}

# --------------------------------------------------------------------------------
# 3. Allometry within each clade
# --------------------------------------------------------------------------------
# A global slope can be produced entirely by between-clade differences in body
# size, with no clade showing that slope internally. Fitting within clades
# separates the two.
if ("clade" %in% colnames(POSTCRANIAL)) {
  d <- POSTCRANIAL[!is.na(POSTCRANIAL$clade), ]

  by_clade <- do.call(rbind, lapply(unique(d$clade), function(cl) {
    s <- d[d$clade == cl, ]
    if (sum(is.finite(s$SkeletalSize)) < MIN_FIT) return(NULL)
    res <- do.call(rbind, lapply(GM_ELEMENTS, fit_element, d = s))
    if (is.null(res)) return(NULL)
    cbind(clade = cl, res)
  }))

  if (!is.null(by_clade)) {
    by_clade$p_adjusted_BH <- signif(p.adjust(by_clade$p_vs_isometry, "BH"), 4)
    cat(sprintf("Within-clade allometry: %d clades with n >= %d\n",
                length(unique(by_clade$clade)), MIN_FIT))
    # Show a fixed number of decimals so every row lines up.
    by_clade$slope         <- format_fixed(by_clade$slope, 4)
    by_clade$ci_lower      <- format_fixed(by_clade$ci_lower, 4)
    by_clade$ci_upper      <- format_fixed(by_clade$ci_upper, 4)
    by_clade$r_squared     <- format_fixed(by_clade$r_squared, 4)
    by_clade$p_adjusted_BH <- format_fixed(by_clade$p_adjusted_BH, 4)
    write_supp(by_clade, "S13_skeletal_allometry_by_clade")
  } else {
    cat("No clade reached n >=", MIN_FIT, "for within-clade allometry.\n")
  }
}

cat("\nNOTE: these regressions treat specimens as independent. With Pagel's\n")
cat("lambda near 0.7, a PGLS is the correct test for any slope reported.\n\n")
