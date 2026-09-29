transform_asc_file <- function(input_path, output_path) {
  ncols <- NA
  nrows <- NA
  
  # Safety wrapper similar to Python try/except
  tryCatch({
    infile  <- file(input_path, open = "r")
    outfile <- file(output_path, open = "w")
    
    on.exit({
      close(infile)
      close(outfile)
    }, add = TRUE)
    
    # Read and parse the first 6 header lines
    for (i in 1:6) {
      line <- trimws(readLines(infile, n = 1))
      parts <- strsplit(line, "\\s+")[[1]]
      
      if (i == 1) {        # ncols line
        ncols <- parts[2]
      } else if (i == 2) { # nrows line
        nrows <- parts[2]
      }
    }
    
    if (is.na(ncols) || is.na(nrows)) {
      stop("Could not find ncols or nrows in the ASC header")
    }
    
    # Write simplified header
    writeLines(paste(ncols, nrows), outfile)
    
    # Copy the remainder of the file (raster data)
    repeat {
      line <- readLines(infile, n = 1)
      if (length(line) == 0) break
      writeLines(line, outfile)
    }
    
    return(output_path)
    
  }, error = function(e) {
    message("Error transforming ASC file: ", e$message)
    
    if (file.exists(output_path)) {
      try(unlink(output_path), silent = TRUE)
    }
    
    return(NULL)
  })
}


process_folder <- function(asc_folder, output_folder) {
  # Check source folder
  if (!dir.exists(asc_folder)) {
    message("Error: Folder '", asc_folder, "' does not exist")
    return(invisible(NULL))
  }
  
  # Create output folder if needed
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Recursively find .asc files
  asc_files <- list.files(
    asc_folder,
    pattern = "\\.asc$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(asc_files) == 0) {
    message("No .asc files found in '", asc_folder, "' or its subfolders")
    return(invisible(NULL))
  }
  
  message("Found ", length(asc_files), " .asc files to process")
  
  for (asc_file in asc_files) {
    # Relative path
    rel_path <- substring(asc_file, nchar(asc_folder) + 2)
    
    rel_dir  <- dirname(rel_path)
    filename <- basename(rel_path)
    
    # Output directory
    output_dir <- if (rel_dir == ".") {
      output_folder
    } else {
      file.path(output_folder, rel_dir)
    }
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Output filename
    base_name <- tools::file_path_sans_ext(filename)
    output_file <- file.path(output_dir, paste0(base_name, ".txt"))
    
    # Transform
    transform_asc_file(asc_file, output_file)
  }
}


save_raster_as_input_format <- function(raster, output_file) {
  temp_file <- tempfile(pattern = "map_", fileext = ".asc")
  writeRaster(raster, temp_file,
              datatype = "INT2S", overwrite = TRUE, NAflag = -9999)
  transform_asc_file(temp_file, output_file)
  unlink(temp_file)
}
