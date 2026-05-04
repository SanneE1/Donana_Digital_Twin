library(terra)
library(tidyverse)
library(tidyterra)
library(MODISTools)
library(lubridate)
library(ows4R)

#----------------------------------------------------------------------------------------------------------------------------------
# NDVI
#----------------------------------------------------------------------------------------------------------------------------------

Dborder <- vect(DT_borders_path)
base_map <- rast("data/GIS_maps/Rabbit_HabitatMap_500_Donana_Fordham_2013.asc")
transects <- vect("data/original_data/Rabbit_donana_KAI_PacoCarro/Transect_oryctolagus.kml")
transects_pr <- project(transects, crs(Dborder))

# Download NDVI for a specific location and time period   - VNP13A1 = VIIRS 500m
print("downloading NDVI raster, will take a few minutes")

DT_84 <- project(Dborder, crs(transects))
DT_ext <- ext(DT_84)

ndvi_data <- mt_subset(
  product = "MOD13Q1",  # This is the 16-day 500m product
  lon = mean(DT_ext[c(1,2)]),           # Your latitude
  lat = mean(DT_ext[c(3,4)]),         # Your longitude
  band = "250m_16_days_NDVI",
  start = "2004-01-01",
  km_lr = 20,            # km left and right
  km_ab = 25,            # km above and below
  site_name = "DT"
)

ndvi_raster <- mt_to_terra(ndvi_data, reproject = TRUE)
writeRaster(ndvi_raster, filename = "data/GIS_maps/ndvi.tif", overwrite = TRUE)
