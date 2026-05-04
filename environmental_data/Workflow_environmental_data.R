library(dplyr)
library(terra)

# Set paths---------------------------------------------------------------------
data_dir = file.path("environmental_data", "data")
border_path = file.path(data_dir, "Donana_DT_border", "Limite_Don╠âana.shp")
template_path = file.path(data_dir, "template_raster_500.tif")
dem_path = file.path(data_dir, "Copernicus_GLO90_Europe_250m.tif")

# load functions----------------------------------------------------------------
# this file just assigns two objects: API_User and API_Key but stored in data to keep it private
source(file.path(data_dir, "CDS_API_info.R"))
source(file.path("environmental_data", "functions", "Download_and_Krig_temperature.R"))


#-------------------------------------------------------------------------------
# Create the raster template
#-------------------------------------------------------------------------------
border_vec <- vect(border_path)
border_vec <- buffer(border_vec, 1000)

# Create template rasters with 500m resolution
DT_template <- rast(
  xmin = xmin(border_vec), 
  xmax = xmax(border_vec),
  ymin = ymin(border_vec), 
  ymax = ymax(border_vec),
  resolution = c(500, 500),
  crs = crs(border_vec)
)
values(DT_template) <- 0

writeRaster(DT_template, template_path, 
            overwrite = TRUE, NAflag = -9999)

#-------------------------------------------------------------------------------
# Download and Krig temperature data
#-------------------------------------------------------------------------------

download_and_krig_temp(template_path = template_path, dem_path = dem_path, 
                       output_dir = file.path(data_dir, "CDS"))


#-------------------------------------------------------------------------------
# Download precipitation data
#-------------------------------------------------------------------------------

cat('still need to figure out how to do this only for the months not yet downloaded')

if(!file.exists(file.path(data_dir, "CDS", "precipitation.grib"))){
  cat('downloading precipitation data - this might take a while')
  # use_python("C:/Users/z1512834z/AppData/Local/Programs/Python/Python314/python.exe")
  source_python('python/download_precipitation_CDS.py')
}

#-------------------------------------------------------------------------------
# Download and clean NDVI
#-------------------------------------------------------------------------------

# More info in ndvi/python/NDVI-Download-Pipeline/README.md

# I had problems with a conflicting proj db. Creating a clean pyton environment first to avoid conflict
# env_dir <- "python_env"
# python_exe <- file.path(env_dir, "Scripts", "python.exe")
# 
# if (!file.exists(python_exe)) {
#   system(sprintf('python -m venv "%s"', env_dir))
#   
#   system(sprintf('"%s" -m pip install --upgrade pip', python_exe))
#   system(sprintf('"%s" -m pip install rasterio pyproj', python_exe))
#   
#   # Installing required modules
#   system2(python_exe, args = c("-m", "pip", "install", "-r", "environmental_data/functions/NDVIDownloadPipeline/requirements.txt"))
#   
# }
# 
# proj_path <- system2(
#   python_exe,
#   args = c("-c", shQuote("import pyproj; print(pyproj.datadir.get_data_dir())")),
#   stdout = TRUE
# )
# proj_path <- trimws(proj_path)
# 
# system2(
#   python_exe,
#   args = c("environmental_data/functions/NDVIDownloadPipeline/pipeline.py", 
#            "--output-dir", "environmental_data/data/NDVI/",
#            "--shape-path", "environmental_data/data/Donana_DT_border/Limite_Don╠âana.shp")
# )






