parse_nts <- function(filepath) {
  # Read the outline file and trim extra spaces from each line.
  lines <- readLines(filepath, warn = FALSE)
  lines <- trimws(lines)

  # Read the first line to get the basic file layout.
  header  <- as.integer(strsplit(lines[1], "\\s+")[[1]])
  npts    <- header[1]   # points per outline (100)
  ndim    <- header[2]   # dimensions (2)
  nspec   <- header[3]   # number of specimens (197)

  # Read the specimen names listed near the top of the file.
  names_raw <- lines[3:(2 + nspec)]

  # Find where each specimen's coordinate block starts.
  block_starts <- grep("^'", lines)

  if (length(block_starts) != nspec) {
    warning(sprintf("Expected %d data blocks, found %d", nspec, length(block_starts)))
  }

  outlines <- vector("list", length(block_starts))

  for (i in seq_along(block_starts)) {
    # Take the specimen name from the block label.
    blk_name <- sub("^'", "", lines[block_starts[i]])

    # Read the x and y coordinates for this outline.
    coord_start <- block_starts[i] + 1
    coord_end   <- coord_start + npts - 1
    coord_lines <- lines[coord_start:coord_end]

    coords <- matrix(NA_real_, nrow = npts, ncol = 2)
    for (j in seq_len(npts)) {
      vals <- as.numeric(strsplit(coord_lines[j], "\\s+")[[1]])
      if (length(vals) >= 2) coords[j, ] <- vals[1:2]
    }

    outlines[[i]] <- data.frame(x = coords[, 1], y = coords[, 2])
    names(outlines)[i] <- blk_name
  }

  cat(sprintf("  Parsed: %d specimens × %d points\n", length(outlines), npts))
  return(list(outlines = outlines, npts = npts, nspec = nspec,
              specimen_names = names_raw))
}

nts_data    <- parse_nts("data/Pterosaur_Outlines_Clean.nts")

outlines_list <- nts_data$outlines
