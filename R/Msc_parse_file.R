parse_nts <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- trimws(lines)
  
  # --- Header ----------------------------------------------------------------
  header  <- as.integer(strsplit(lines[1], "\\s+")[[1]])
  npts    <- header[1]   # points per outline (100)
  ndim    <- header[2]   # dimensions (2)
  nspec   <- header[3]   # number of specimens (197)
  
  # --- Names (lines 3 to 2+nspec, line 2 = "names") -------------------------
  names_raw <- lines[3:(2 + nspec)]
  
  # --- Data blocks -----------------------------------------------------------
  # Lines starting with ' mark the beginning of each specimen block
  block_starts <- grep("^'", lines)
  
  if (length(block_starts) != nspec) {
    warning(sprintf("Expected %d data blocks, found %d", nspec, length(block_starts)))
  }
  
  outlines <- vector("list", length(block_starts))
  
  for (i in seq_along(block_starts)) {
    # Specimen name from block header (remove leading ')
    blk_name <- sub("^'", "", lines[block_starts[i]])
    
    # Coordinate lines
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


nts_data    <- parse_nts("Pterosaur_Outlines_Clean.nts")

outlines_list <- nts_data$outlines