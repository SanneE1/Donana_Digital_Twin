

create_detectionP_map <- function(Cov_list, 
                                  obs_p = obs_p, 
                                  obs_p_ndvi = obs_p_ndvi,
                                  output_dir, start_year, end_year) {
  
  ndvi_files  <- list.files(Cov_list$ndvi_map,  full.names = TRUE)
  ndvi_all <- rast(ndvi_files)
  ndvi_time <- floor_date(time(ndvi_all), "month")
  

  for (year in start_year:end_year) {
    for (month in 1:12){
      
      ndvi <- subset(ndvi_all, ndvi_time == as.Date(paste(year, sprintf("%02d", month), "01", sep = "-")))
      
      DPrast <- terra::app(ndvi, function(x) {boot::inv.logit(obs_p + obs_p_ndvi * x)})
      
      writeRaster(DPrast, file.path(output_dir, paste0(year, "_", sprintf("%02d", month), ".asc")), 
                  datatype = "FLT4S", overwrite = TRUE, NAflag = -9999)
      
      rm(DPrast)
      gc()
    }
  }
}
