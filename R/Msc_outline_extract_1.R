extract_simple_outlines_fuzzy <- function(image_folder, 
                                          specimen_names, 
                                          output_file = "Pterosaur_Outlines",
                                          threshold = 50,
                                          num_points = 100) {
  
  image_files <- list.files(image_folder, 
                            pattern = "\\.(jpg|jpeg|png)$", 
                            ignore.case = TRUE, 
                            full.names = FALSE)
  
  image_paths <- list.files(image_folder, 
                            pattern = "\\.(jpg|jpeg|png)$", 
                            ignore.case = TRUE, 
                            full.names = TRUE)
  
  cat(paste("Found", length(image_files), "images in", image_folder, "\n"))
  cat(paste("Processing", length(specimen_names), "specimens\n\n"))
  
  find_matching_image <- function(specimen_name, image_files) {
    specimen_clean <- gsub(" ", "_", specimen_name)
    specimen_clean <- gsub("\\s+", "_", specimen_clean)
    specimen_clean <- trimws(specimen_clean)
    
    exact_match <- which(grepl(specimen_clean, image_files, ignore.case = TRUE))
    if (length(exact_match) > 0) return(exact_match[1])
    
    species_parts <- strsplit(specimen_name, " ")[[1]]
    if (length(species_parts) >= 2) {
      genus <- trimws(species_parts[1])
      species <- trimws(species_parts[2])
      
      patterns <- c(
        paste0(genus, "_", species),
        paste0("Out_", genus, "_", species),
        paste0(genus, " ", species),
        genus
      )
      
      for (pattern in patterns) {
        match <- which(grepl(pattern, image_files, ignore.case = TRUE))
        if (length(match) > 0) return(match[1])
      }
    }
    
    return(NA)
  }
  
  all_outlines <- list()
  success_log <- data.frame(
    ID = character(),
    Filename = character(),
    Success = integer(),
    NumPoints = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(specimen_names)) {
    specimen_name <- specimen_names[i]
    image_index <- find_matching_image(specimen_name, image_files)
    
    if (is.na(image_index)) {
      if (i %% 20 == 1) {
        cat(paste("⚠️  [", i, "/", length(specimen_names), "] Image not found for:", specimen_name, "\n"))
      }
      success_log <- rbind(success_log, data.frame(
        ID = specimen_name, Filename = NA, Success = 0, NumPoints = 0
      ))
      next
    }
    
    image_file <- image_paths[image_index]
    
    if (i %% 10 == 0 || i <= 10) {
      cat(paste("Processing [", i, "/", length(specimen_names), "]:", specimen_name, "... "))
    }
    
    tryCatch({
      img <- image_read(image_file)
      img_gray <- image_convert(img, colorspace = "gray")
      img_matrix <- as.integer(image_data(img_gray, channels = "gray"))
      binary_mask <- img_matrix > threshold
      
      coords <- which(binary_mask, arr.ind = TRUE)
      
      if (nrow(coords) < 10) {
        if (i %% 10 == 0 || i <= 10) cat("Failed (too few points)\n")
        success_log <- rbind(success_log, data.frame(
          ID = specimen_name, Filename = basename(image_file), Success = 0, NumPoints = nrow(coords)
        ))
        next
      }
      
      hull_indices <- chull(coords[, 2], coords[, 1])
      boundary <- coords[hull_indices, ]
      outline_df <- data.frame(x = boundary[, 2], y = boundary[, 1])
      
      if (nrow(outline_df) < 3) {
        if (i %% 10 == 0 || i <= 10) cat("Failed (outline too small)\n")
        success_log <- rbind(success_log, data.frame(
          ID = specimen_name, Filename = basename(image_file), Success = 0, NumPoints = nrow(outline_df)
        ))
        next
      }
      
      distances <- c(0, cumsum(sqrt(diff(outline_df$x)^2 + diff(outline_df$y)^2)))
      total_distance <- max(distances)
      target_distances <- seq(0, total_distance, length.out = num_points)
      
      resampled_x <- approx(distances, outline_df$x, xout = target_distances)$y
      resampled_y <- approx(distances, outline_df$y, xout = target_distances)$y
      
      outline_resampled <- data.frame(x = resampled_x, y = resampled_y)
      all_outlines[[specimen_name]] <- outline_resampled
      
      if (i %% 10 == 0 || i <= 10) {
        cat(paste("Success (", nrow(outline_resampled), "points)\n"))
      }
      
      success_log <- rbind(success_log, data.frame(
        ID = specimen_name, 
        Filename = basename(image_file), 
        Success = 1, 
        NumPoints = nrow(outline_resampled)
      ))
      
    }, error = function(e) {
      if (i %% 10 == 0 || i <= 10) {
        cat(paste("Error:", e$message, "\n"))
      }
      success_log <- rbind(success_log, data.frame(
        ID = specimen_name, Filename = basename(image_file), Success = 0, NumPoints = 0
      ))
    })
  }
  
  if (length(all_outlines) > 0) {
    num_specimens <- length(all_outlines)
    num_points_check <- nrow(all_outlines[[1]])
    
    cat("\nSaving outlines to NTS format...\n")
    
    con <- file(paste0(output_file, ".nts"), "w")
    
    writeLines(c(
      paste(num_points_check, 2, num_specimens),
      "names",
      paste(names(all_outlines), collapse = "\n")
    ), con)
    
    for (specimen_name in names(all_outlines)) {
      coords <- all_outlines[[specimen_name]]
      writeLines(paste0("'", specimen_name), con)
      
      for (j in 1:nrow(coords)) {
        writeLines(paste(coords$x[j], coords$y[j]), con)
      }
    }
    
    close(con)
    cat(paste("✓ Saved", length(all_outlines), "outlines to", paste0(output_file, ".nts"), "\n"))
  }
  
  write.table(success_log, 
              paste0(output_file, "_Success.txt"), 
              row.names = FALSE, 
              sep = "\t", 
              quote = FALSE)
  return(list(outlines = all_outlines, success_log = success_log))
}

# Run extraction
# Raw silhouette images are not tracked in this repository (large binary assets).
# Place them locally in image_folder to re-run extraction from scratch; otherwise
# the pre-extracted outlines already shipped in data/Pterosaur_Outlines_Clean.nts
# are used by Msc_outline_loading.R and this step is skipped.
image_folder <- "./data/pteros_out_finale/"
specimen_names <- complete_data$SPECIES..Pegas.

if (dir.exists(image_folder)) {
  extraction_results <- extract_simple_outlines_fuzzy(
    image_folder = image_folder,
    specimen_names = specimen_names,
    output_file = "data/Pterosaur_Outlines_Clean",
    threshold = 50,
    num_points = 100
  )
} else {
  cat("Image folder '", image_folder, "' not found - skipping outline extraction ",
      "and using the pre-extracted data/Pterosaur_Outlines_Clean.nts instead.\n", sep = "")
}