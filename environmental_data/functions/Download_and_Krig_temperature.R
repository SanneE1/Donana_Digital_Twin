# devtools::install_github("https://github.com/ErikKusch/KrigR", ref = "Development")

library(KrigR)
library(terra)
library(lubridate)

download_and_krig_temp <- function(template_path, dem_path, output_dir){
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = T)

  template_rast <- rast(template_path)
  hab_rast <- project(template_rast,"epsg:4326")
  
  dem_rast <- rast(dem_path)
  dem_rast <- crop(dem_rast, hab_rast)
  
  cat('Starting Temperature downloading and downscaling')
  
  start <- ymd("2005-01-01")
  end   <- floor_date(Sys.Date(), "month") - days(1)
  
  months_seq <- seq(start, end, by = "month")
  
  for (m in months_seq) {
    
    start_date <- format(m, "%Y-%m-01 00:00")
    end_date   <- format(ceiling_date(m, "month") - days(1), "%Y-%m-%d 23:00")
    
    if(!file.exists(file.path(output_dir, paste0("CDSTemp_Raw_", m, ".nc")))){
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
        Dir = output_dir, # where to store the data
        FileName = paste0("CDSTemp_Raw_", m, ".nc"), # what to call the resulting file
        ## API User Credentials
        API_User = API_User,
        API_Key = API_Key
      )} else {
        cat('reading existing nc file for year:', year, 'half:', half, fill = T)
        CDSTemp_Raw <- rast(file.path(output_dir, paste0("CDSTemp_Raw_", m, ".nc")))
      }
    
    CDSTemp_Raw <- project(CDSTemp_Raw, hab_rast)
    
    if(!exists("Covs_ls")) {
      Covs_ls <- CovariateSetup(
        Training = CDSTemp_Raw,
        Target = res(hab_rast)[1],
        Covariates = dem_rast,
        Source = "Drive",
        Dir = output_dir,
        Keep_Global = TRUE
      )
    } 
    
    if(!file.exists(name)){
      CDSTemp_Krig <- Kriging(
        Data = crop(CDSTemp_Raw[[i]], terra::ext(Covs_ls$Training)), # data we want to krig as a raster object
        Covariates_training = Covs_ls[[1]], # training covariate as a raster object
        Covariates_target = Covs_ls[[2]], # target covariate as a raster object
        # Equation = "output_hh", # the covariate(s) we want to use
        nmax = 40, # degree of localisation
        Cores = 3, # we want to krig using three cores to speed this process up
        FileName = m, # the file name for our full kriging output
        Dir = output_dir # which directory to save our final input in
      )
    } 
  }
}



