# devtools::install_github("https://github.com/ErikKusch/KrigR", ref = "Development")

library(KrigR)
library(terra)
library(lubridate)

# this file just assigns two objects: API_User and API_Key but stored in data to keep it private
source("data/CDS_API_info.R")  

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "data", "original_data", "CDS") # folder path for data

## create directories, if they don't exist yet
if (!dir.exists(Dir.Data)) dir.create(Dir.Data)
if (!dir.exists(file.path(Dir.Data, "downscaled_monthly_files_donana"))) {dir.create(file.path(Dir.Data, "downscaled_monthly_files_donana"), recursive = T)}
if (!dir.exists(file.path(Dir.Data, "monthly_files_donana"))) {dir.create(file.path(Dir.Data, "monthly_files_donana"), recursive = T)}
if (!dir.exists(file.path("data", "model_input", "climate_historic"))) {dir.create(file.path("data", "model_input", "climate_historic"), recursive = T)}


hab_rast <- rast("data/GIS_maps/Rabbit_HabitatMap_500_Donana_Fordham_2013.asc")
hab_rast <- project(hab_rast,"epsg:4326")
hab_rast <- extend(hab_rast, 100)

dem_rast <- rast("data/original_data/Copernicus_GLO90_Europe_250m.tif")
dem_rast <- crop(dem_rast, hab_rast)

coords_donana <- read.table("data/coordinates_DT_EPSG4326.txt")[,c(1,2)]

cat('Starting Temperature downloading and downscaling')

for (year in 2005:2026) {
  for (half in 1:2) {
    if (half == 1) {
      start_date <- paste0(year, "-01-01 00:00")
      end_date   <- paste0(year, "-06-30 23:00")
    } else {
      start_date <- paste0(year, "-07-01 00:00")
      end_date   <- paste0(year, "-12-31 23:00")
    }
    
    if(!file.exists(file.path(Dir.Data, paste0("CDSTemp_Raw_", year, "_", half, ".nc")))){
      CDSTemp_Raw <- CDownloadS(
        ## Variable and Data Product
        Variable = "2m_temperature", # this is air temperature
        DataSet = "reanalysis-era5-land", # data product from which we want to download
        ## Time-Window
        DateStart = start_date, # date at which time window opens
        DateStop = end_date, # date at which time window terminates
        TZone = "CET", # European Central Time to align with our study region
        ## Temporal Aggregation
        TResolution = "month", # we want daily aggregates
        TStep = 1, # we want aggregates of 1 day each
        ## Spatial Limiting
        Extent = hab_rast, # our rectangular bounding box
        ## File Storing
        Dir = Dir.Data, # where to store the data
        FileName = paste("CDSTemp_Raw", year, half, sep = "_"), # what to call the resulting file
        ## API User Credentials
        API_User = API_User,
        API_Key = API_Key
      )} else {
        cat('reading existing nc file for year:', year, 'half:', half, fill = T)
        CDSTemp_Raw <- rast(file.path(Dir.Data, paste0("CDSTemp_Raw_", year, "_", half, ".nc")))
      }
    
    CDSTemp_Raw <- project(CDSTemp_Raw, hab_rast)
    
    if(!exists("Covs_ls")) {
      Covs_ls <- CovariateSetup(
        Training = CDSTemp_Raw,
        Target = res(hab_rast)[1],
        Covariates = dem_rast,
        Source = "Drive",
        Dir = Dir.Data,
        Keep_Global = TRUE
      )
    } 
    
    for(i in c(1:length(time(CDSTemp_Raw)))){
      name <- paste(lubridate::year(time(CDSTemp_Raw[[i]])), lubridate::month(time(CDSTemp_Raw[[i]])), "Krigged.nc", sep = "_")
      
      if(!file.exists(name)){
        CDSTemp_Krig <- Kriging(
          Data = crop(CDSTemp_Raw[[i]], terra::ext(Covs_ls$Training)), # data we want to krig as a raster object
          Covariates_training = Covs_ls[[1]], # training covariate as a raster object
          Covariates_target = Covs_ls[[2]], # target covariate as a raster object
          # Equation = "output_hh", # the covariate(s) we want to use
          nmax = 40, # degree of localisation
          Cores = 3, # we want to krig using three cores to speed this process up
          FileName = name, # the file name for our full kriging output
          Dir = Dir.Data # which directory to save our final input in
        )
      } else {
        CDSTemp_Krig <- rast(name)
      }
      
      if(file.exists(file.path(Dir.Data, "downscaled_monthly_files_donana",
                               paste0("tas_", year(time(CDSTemp_Raw[[i]])), "_", sprintf("%02d", month(time(CDSTemp_Raw[[i]]))), ".txt")))){
        next 
      } else {
        temp_donana <- extract(CDSTemp_Krig$Prediction, coords_donana)[,2] - 273.15
        
        write.table(temp_donana,
                    file.path(Dir.Data, "downscaled_monthly_files_donana",
                              paste0("tas_", year(time(CDSTemp_Raw[[i]])), "_", sprintf("%02d", month(time(CDSTemp_Raw[[i]]))), ".txt")),
                    quote = FALSE, row.names = FALSE, col.names = F)
        
      }
      
      
      
    }
    
  }
}


cat('Starting Precipitation downloading')

if(!file.exists(file.path("data", "original_data", "CDS", "precipitation.grib"))){

  cat('downloading precipitation data - this might take a while')
  library(reticulate)
  use_python("C:/Users/z1512834z/AppData/Local/Programs/Python/Python314/python.exe")
  source_python('python/download_precipitation_CDS.py')

}

precip <- rast(file.path("data", "original_data", "CDS", "precipitation.grib"))


for(i in c(1:length(time(precip)))){

  date = time(precip[[i]])
  Precip_donana <- extract(precip[[i]], coords_donana)
  # Convert average daily precipitation in M to total monthly precipication in mm
  Precip_donana <- Precip_donana * 1000 * as.integer(lubridate::days_in_month(date))

  write.table(Precip_donana,
              file.path("data", "original_data", "CDS", "monthly_files_donana",
                        paste0("pr_", year(date), "_", sprintf("%02d", month(date)), ".txt")),
              quote = FALSE, row.names = FALSE, col.names = F)

}

# ---------------------------------------------------------------------------------------------
# Format to rabbit variables
# ---------------------------------------------------------------------------------------------

source("R/Function_format_climate_data_for_rabbit.R")

# ERA5-Land 
calculate_BM_cDM(tas_files = file.path(Dir.Data, "downscaled_monthly_files_donana"),
                 pr_files = file.path(Dir.Data, "monthly_files_donana"),
                 year_min = 2006,
                 year_max = 2024,
                 result_dir = "data/model_input/climate_historic_COP/",
                 coord_file = "data/coordinates_DT_EPSG4326.txt",
                 temperature_type = "ERA5")

