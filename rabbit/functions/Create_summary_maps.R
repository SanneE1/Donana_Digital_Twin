calculate_matrix_statistics <- function(folder_paths, template_raster, confidence_level = 0.95) {

  # Get the map dirs
  map_dirs <- list.dirs(folder_paths, recursive = T)
  map_dirs <- map_dirs[grepl("maps", map_dirs) ]
  
  # Get list of CSV files from the first folder
  first_folder <- file.path(map_dirs[1])
  csv_files <- sort(list.files(first_folder, pattern = "\\.csv$", full.names = TRUE))
  
  if (length(csv_files) == 0) {
    cat(sprintf("No CSV files found in %s\n", first_folder))
    return(invisible(NULL))
  }
  
  cat(sprintf("Found %d CSV files\n", length(csv_files)))
  cat(sprintf("Processing %d folders\n", length(map_dirs)))
  
  # Output: named list, one entry per CSV file
  results <- list()
  
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    base_name <- sub("\\.csv$", "", filename)
    
    # Collect all matrices for this file
    matrices <- list()
    
    for (folder in map_dirs) {
      file_path <- file.path(folder, filename)
      if (file.exists(file_path)) {
        matrices[[length(matrices) + 1]] <- as.matrix(read.csv(file_path, header = FALSE))
      } else {
        cat(sprintf("  Warning: %s not found, skipping\n", file_path))
      }
    }
    
    if (length(matrices) == 0) {
      cat(sprintf("  No valid matrices found for %s, skipping\n", filename))
      next
    }
    
    # Stack into 3D array (rows x cols x n_folders)
    n_samples <- length(matrices)
    n_rows    <- nrow(matrices[[1]])
    n_cols    <- ncol(matrices[[1]])
    matrices_array <- array(unlist(matrices), dim = c(n_rows, n_cols, n_samples))
    
    # Compute per-cell statistics
    mean_matrix  <- apply(matrices_array, c(1, 2), mean)
    std_matrix   <- apply(matrices_array, c(1, 2), sd)
    se_matrix    <- std_matrix / sqrt(n_samples)
    t_value      <- qt((1 + confidence_level) / 2, df = n_samples - 1)
    margin_error <- t_value * se_matrix
    lower_ci     <- mean_matrix - margin_error
    upper_ci     <- mean_matrix + margin_error
    
    # Helper: convert a matrix to a SpatRaster using the template
    to_raster <- function(mat) {
      r <- template_raster
      values(r) <- as.vector(t(mat))  # t() to match row-major CSV order
      r
    }
    
    results[[base_name]] <- list(
      mean     = to_raster(mean_matrix),
      std      = to_raster(std_matrix),
      lower_ci = to_raster(lower_ci),
      upper_ci = to_raster(upper_ci)
    )
  }
  
  return(results)
}
