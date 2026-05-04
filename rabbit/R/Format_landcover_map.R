# Load required libraries
library(terra)
library(tidyverse)
library(tidyterra)
library(MODISTools)
library(lubridate)
library(ows4R)

corine_raster_path = "data/original_data/U2018_CLC2018_V2020_20u1.tif"
DT_borders_path = "data/original_data/Donana_DT_border/Limite_Don╠âana.shp"
output_folder = "data/GIS_maps/"

# Load DT boundaries
DT_vect <- vect(DT_borders_path)
DT_vect_B <- buffer(DT_vect, 1000)

#----------------------------------------------------------------------------------------------------------------------------------
# CORINE 
#----------------------------------------------------------------------------------------------------------------------------------

corine_raster <- rast(corine_raster_path)

corine_raster <- crop(corine_raster, ext(c(2546600, 3075620, 1480300, 1921470))) # shortcut to shorten the project

corine_raster <- project(corine_raster, crs(DT_vect))
corine_raster <- crop(corine_raster, DT_vect_B)

# Reclassification table --- Based on Fordham 2013 Table S3
reclass_mat <- as.matrix(data.frame(
  old = c(1:44, 48),
  new = c(rep(0,9), 15, 0,15,0,0, 15,15,15,30,15,15,15,15,0,0,15,30,15,15,30,15,0,30,15,rep(0,12))
))

reclas_donana <- classify(corine_raster, reclass_mat)

# Resize to 500x500m WITHOUT using project()
print("Resizing to 500x500m resolution...")

# Create template rasters with 500m resolution
DT_template <- rast(
  xmin = xmin(reclas_donana), 
  xmax = xmax(reclas_donana),
  ymin = ymin(reclas_donana), 
  ymax = ymax(reclas_donana),
  resolution = c(500, 500),
  crs = crs(reclas_donana)
)

# Resample instead of project
DT_rast <- resample(reclas_donana, DT_template, method = "modal")
DT_rast <- mask(DT_rast, DT_vect_B)

# Set the output data type to Int16
writeRaster(DT_rast, file.path(output_folder, "Rabbit_HabitatMap_500_Donana_Fordham_2013.asc"), 
            datatype = "INT2S", overwrite = TRUE, NAflag = -9999)

rm(list = ls())



ndvi_raster <- project(ndvi_raster, crs(base_map))
ndvi_raster <- resample(ndvi_raster, base_map, method = "mean")

dates <- as.Date(names(ndvi_raster))

yrs <- unique(year(dates))
mnts <- c(1:12)

for(y in yrs){
  for(m in mnts){
    try({
      t_rast <- subset(ndvi_raster, which(month(dates) == m & year(dates) == y))
      
      t_mean <- app(t_rast, mean, na.rm = T)
      # plot(t_mean)
      writeRaster(t_mean, 
                  file.path("data", "GIS_maps", "ndvi_monthly_maps", paste0(y, "_", sprintf("%02d", m), ".asc")),
                  overwrite = T)
    })
  }
}


#----------------------------------------------------------------------------------------------------------------------------------
# Flooding maps
#----------------------------------------------------------------------------------------------------------------------------------

fld_files <- list.files(file.path("data", "original_data", "LAST_maps", "Inundacion"), full.names = T)

dates <- stringr::str_extract(fld_files, "\\d{8}")
## Below round the date of the flood to the nearest month, so instead of taking floods that happen throughout the month, the flood is taken into
## account at the load of the closest "new month"
dates <- as.Date(dates, format = "%Y%m%d")
yrs <- unique(year(dates))
yrs <- yrs[which(yrs > 1999)]
mnts <- c(1:12)

for(y in yrs){
  for(m in mnts){
    try({
      
      r_files <- fld_files[which(year(dates) == y & month(dates) == m)]
      
      if (length(r_files) == 0) next
      
      r_list <- lapply(r_files, function(f) {
        r <- rast(f)
        r <- project(r, crs(base_map))
        r <- resample(r, base_map, method = "max")
        return(r)
      })
      
      r_max <- round(max(rast(r_list), na.rm = TRUE), digits = 0) 
      
      writeRaster(r_max, 
                  file.path("data", "GIS_maps", "LAST_floodmaps", paste0(y, "_", sprintf("%02d", m), ".asc")),
                  overwrite = T, NAflag = 0)
    })
  }
}










