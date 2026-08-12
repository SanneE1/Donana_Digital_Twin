library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(reticulate)
library(data.table)  
library(lubridate)
library(xgboost)     
library(zoo)         

# Next step (version?) would be to use
# precipitation and temperature stacks 
# instead of the warm/dry options

ndvi_forecast_function <- function(info_dir = file.path("ndvi", "results"),
                                   temp_choices = rep("warm", 24),
                                   precip_choices = rep("mean", 24)
)
{

  sapply(list.files(file.path("ndvi", "functions"), full.names = T), source)
  
  #-------------------------------------------------------------------------------
  # FORECAST
  #-------------------------------------------------------------------------------
  
  temp_warm <- load_stack_with_time(tif_path = here::here(info_dir, "temp_warm_values.tif"),
                                    time_path = here::here(info_dir, "temp_warm_time.rds"))
  temp_mean <- load_stack_with_time(tif_path = here::here(info_dir, "temp_mean_values.tif"),
                                    time_path = here::here(info_dir, "temp_mean_time.rds"))
  temp_cold <- load_stack_with_time(tif_path = here::here(info_dir, "temp_cold_values.tif"),
                                    time_path = here::here(info_dir, "temp_cold_time.rds"))
  
  precip_wet <- load_stack_with_time(tif_path = here::here(info_dir, "precip_wet_values.tif"),
                                     time_path = here::here(info_dir, "precip_wet_time.rds"))
  precip_mean <- load_stack_with_time(tif_path = here::here(info_dir, "precip_mean_values.tif"),
                                      time_path = here::here(info_dir, "precip_mean_time.rds"))
  precip_dry <- load_stack_with_time(tif_path = here::here(info_dir, "precip_dry_values.tif"),
                                     time_path = here::here(info_dir, "precip_dry_time.rds"))
  
  
  clim_list <- list("temp" = list("warm" = temp_warm, "mean" = temp_mean, "cold" = temp_cold),
                    "precip" = list("wet" = precip_wet, "mean" = precip_mean, "dry" = precip_dry))
  
  ndvi_last_date <- local({
    r <- load_stack_with_time(
      tif_path = here::here(info_dir, "ndvi_stack.tif"),
      time_path = here::here(info_dir, "ndvi_time.rds")
    )
    max(time(r))
  })
  
  ndvi_forecast <- shiny_forecast_function(temp_choices = temp_choices,
                                           precip_choices = precip_choices, 
                                           future_climate = clim_list, 
                                           ndvi_last_date = ndvi_last_date,
                                           data_dir = info_dir 
  )
  
  return(ndvi_forecast)
}

ndvi_forecast_function(temp_choices = rep("warm", 24), 
                       precip_choices = rep("average", 24))
ndvi_forecast_function(temp_choices = rep("cold", 24), 
                       precip_choices = rep("average", 24))
ndvi_forecast_function(temp_choices = rep("average", 24), 
                       precip_choices = rep("wet", 24))
ndvi_forecast_function(temp_choices = rep("average", 24), 
                       precip_choices = rep("dry", 24))

ndvi_forecast_function(temp_choices = rep("cold", 24), 
                       precip_choices = rep("wet", 24))
ndvi_forecast_function(temp_choices = rep("warm", 24), 
                       precip_choices = rep("dry", 24))
