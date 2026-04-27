
library(dplyr)
library(terra)

csvToRaster <- function(fileName, habitat_raster) {
  
  mat_status <- as.matrix(read.csv(fileName, header = F))
  
  df <- terra::rasterize(mat_status, habitat_raster)
  values(df) <- mat_status
 
  return(df)  
  }

inputMapToRaster <- function(fileName, habitat_raster) {
  
  mat <- as.matrix(read.table(fileName, sep = " ", header = F, skip = 1))
  r <- terra::rasterize(mat, habitat_raster)
  values(r) <- mat
  
  return(r)
}


