

create_K_map <- function(Cov_list, output_dir, alpha, beta_param, K, start_year, end_year) {
  
  ndvi_files  <- list.files(Cov_list$ndvi_map,  full.names = TRUE)
  flood_files <- list.files(Cov_list$flood_map, full.names = TRUE)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ndvi_all <- rast(ndvi_files)
  flood_all <- rast(flood_files)
  
  ndvi_time <- floor_date(time(ndvi_all), "month")
  flood_time <- floor_date(time(flood_all), "month")
  
  for (year in start_year:end_year) {
    for (month in 1:12){
      
      ndvi <- subset(ndvi_all, ndvi_time == as.Date(paste(year, sprintf("%02d", month), "01", sep = "-")))
      flood <- tryCatch(
        subset(flood_all, flood_time == as.Date(paste(year, sprintf("%02d", month), "01", sep = "-"))),
        error = function(e) NULL
      )
      
      Krast <- ifel(ndvi <= 0, 0, ndvi)
      Krast <- terra::app(Krast, function(x) {dbeta(x, alpha, beta_param)})
      
      if(!is.null(flood)){
        Krast <- ifel(flood == 1, 0, Krast)  
      }
      Krast <- round(Krast * K)
      
      save_raster_as_input_format(Krast, file.path(output_dir, 
                                                   paste0(year, "_", sprintf("%02d", month), ".txt")))
      
      }
  }
}
