cat("\n--- POSTCRANIAL DATA PREPARATION ---\n\n")

# Read the limb measurements and turn them into a cleaned size table.
if (!exists("write_supp")) stop("Source Msc_stats_00_setup.R first.")

MAIN_CSV <- "data/pteros_main_data.csv"
if (!file.exists(MAIN_CSV)) stop("'", MAIN_CSV, "' not found in the working directory.")

raw <- read.csv(MAIN_CSV, sep = ";", stringsAsFactors = FALSE)

# Trailing dots come from trailing spaces in the source headers
ELEMENTS <- c(
  Humerus     = "Humerus.",
  Ulna        = "Ulna",
  Pteroid     = "Pteroid",
  McIV        = "Mc.IV",
  WingPh1     = "Wg.ph.1",
  WingPh2     = "Wg.ph.2",
  WingPh3     = "Wg.ph.3",
  WingPh4     = "Wg.ph.4",
  Femur       = "Femur.",
  Tibia       = "Tibia",
  MtIII       = "Mt.III",
  WingLength  = "Wing.length...cm.",
  HindLimb    = "measured.hind.limbs..cm.",
  Torso       ="DSV"       )

missing_el <- ELEMENTS[!ELEMENTS %in% colnames(raw)]
if (length(missing_el) > 0) {
  cat("Columns not found (check for trailing dots):",
      paste(missing_el, collapse = ", "), "\n")
  ELEMENTS <- ELEMENTS[ELEMENTS %in% colnames(raw)]
}
cat("Postcranial elements available:", length(ELEMENTS), "\n")

# --------------------------------------------------------------------------------
# 1. Assemble and clean
# --------------------------------------------------------------------------------
post <- data.frame(species = trimws(raw$SPECIES..Pegas.), stringsAsFactors = FALSE)
for (nm in names(ELEMENTS)) {
  post[[nm]] <- suppressWarnings(as.numeric(raw[[ELEMENTS[[nm]]]]))
}

EL_NAMES <- names(ELEMENTS)

# Zeros are absences, not measurements of length zero: they would send log() to
# -Inf and drag the geometric mean to zero
for (nm in EL_NAMES) post[[nm]][post[[nm]] <= 0] <- NA

cat("\nCompleteness per element:\n")
comp <- data.frame(element = EL_NAMES,
                   n = vapply(EL_NAMES, function(nm) sum(is.finite(post[[nm]])), integer(1)),
                   min = round(vapply(EL_NAMES, function(nm)
                     suppressWarnings(min(post[[nm]], na.rm = TRUE)), numeric(1)), 2),
                   max = round(vapply(EL_NAMES, function(nm)
                     suppressWarnings(max(post[[nm]], na.rm = TRUE)), numeric(1)), 2),
                   stringsAsFactors = FALSE)
print(comp, row.names = FALSE)

# --------------------------------------------------------------------------------
# 2. Implausible values
# --------------------------------------------------------------------------------
# A wrist or limb element longer than the wing itself is a data-entry error, not
# anatomy. Pteroid reaches 504 cm in this database while Mc IV peaks at 125 cm.
flag <- do.call(rbind, lapply(setdiff(EL_NAMES, c("WingLength", "HindLimb")), function(nm) {
  bad <- which(is.finite(post[[nm]]) & is.finite(post$WingLength) &
                 post[[nm]] > post$WingLength)
  if (length(bad) == 0) return(NULL)
  data.frame(species = post$species[bad], element = nm,
             value = post[[nm]][bad], wing_length = post$WingLength[bad],
             stringsAsFactors = FALSE)
}))

if (!is.null(flag)) {
  cat(sprintf("\nIMPLAUSIBLE: %d element measurements exceed the wing length\n",
              nrow(flag)))
  print(head(flag, 15), row.names = FALSE)
  cat("These are set to NA. Check the source rows before publishing.\n")
  write_supp(flag, "S11_implausible_postcranial_values")
  for (i in seq_len(nrow(flag))) {
    post[[flag$element[i]]][post$species == flag$species[i]] <- NA
  }
}

# --------------------------------------------------------------------------------
# 3. Geometric-mean size proxy
# --------------------------------------------------------------------------------
# Use the measured bones to build one overall body-size proxy per specimen.
# Computed on the elements only, excluding the two summary lengths, which are
# sums of the others and would be counted twice.
GM_ELEMENTS <- setdiff(EL_NAMES, c("WingLength", "HindLimb"))

post$SkeletalSize <- exp(rowMeans(log(post[, GM_ELEMENTS]), na.rm = TRUE))
post$n_elements   <- rowSums(is.finite(as.matrix(post[, GM_ELEMENTS])))

# A geometric mean over three bones is not comparable to one over ten
post$SkeletalSize[post$n_elements < 6] <- NA

cat(sprintf("\nSkeletal size proxy computed for %d specimens (>= 6 elements)\n",
            sum(is.finite(post$SkeletalSize))))
print(summary(post$SkeletalSize))

# --------------------------------------------------------------------------------
# 4. Join to the performance data
# --------------------------------------------------------------------------------
POSTCRANIAL <- merge(performance_data_clean, post, by = "species", all.x = TRUE)

cat(sprintf("\nJoined: %d of %d specimens have a skeletal size\n",
            sum(is.finite(POSTCRANIAL$SkeletalSize)), nrow(POSTCRANIAL)))

unmatched <- setdiff(performance_data_clean$species, post$species)
if (length(unmatched) > 0) {
  cat("Unmatched species:", length(unmatched), "\n")
  print(head(unmatched, 10))
}

write_supp(comp, "S11b_postcranial_completeness")
cat("\nAvailable: POSTCRANIAL, EL_NAMES, GM_ELEMENTS\n\n")
