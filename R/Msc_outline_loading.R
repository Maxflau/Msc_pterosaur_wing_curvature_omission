if (file.exists("data/Pterosaur_Outlines_Clean.nts")) {
  
  cat("Loading extracted outlines...\n")
  
  lines <- readLines("data/Pterosaur_Outlines_Clean.nts")
  header <- as.numeric(strsplit(lines[1], " ")[[1]])
  num_points <- header[1]
  num_specimens <- header[3]
  
  name_start <- 3
  name_end <- name_start + num_specimens - 1
  specimen_names_nts <- lines[name_start:name_end]
  
  outlines_list <- list()
  current_line <- name_end + 1
  
  for (i in 1:num_specimens) {
    current_line <- current_line + 1
    coords <- matrix(NA, nrow = num_points, ncol = 2)
    
    for (j in 1:num_points) {
      coord_values <- as.numeric(strsplit(lines[current_line], " ")[[1]])
      coords[j, ] <- coord_values
      current_line <- current_line + 1
    }
    
    outlines_list[[specimen_names_nts[i]]] <- data.frame(x = coords[,1], y = coords[,2])
  }
  
  cat(paste("✓ Loaded", length(outlines_list), "wing outlines\n\n"))
  
} else {
  stop("Outline file not found! Run outline extraction first.")
}