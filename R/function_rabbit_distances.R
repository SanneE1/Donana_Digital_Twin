library(dplyr)
library(terra)
library(sf)
library(lubridate)

# sim_output = "Out_dir" 
# hab_file = "data/GIS_maps/Rabbit_HabitatMap_500_Donana_Fordham_2013.asc" 
# transectsD = "data/original_data/Rabbit_donana_KAI_PacoCarro/Transect_oryctolagus.kml"
# obs_data = read.csv("data/original_data/Rabbit_donana_KAI_PacoCarro/KAI_Rabbit_Night_2024_v1.csv")
# obs_p = 0
# obs_p_ndvi = 0

get_sim_data <- function(sim_output, hab_file, transectsD, obs_data, obs_p, obs_p_ndvi) {
  suppressWarnings(suppressMessages({
  obs_data <- obs_data %>% mutate(year = year(as.Date(Fecha, "%d/%m/%Y")),
                                  month = month(as.Date(Fecha, "%d/%m/%Y")))
  dates <- unique(data.frame(year = obs_data$year, month = obs_data$month))
  
  transects <- vect(transectsD)
  vect_length <- terra::perim(transects)/1000  # Perim gives m, but need km for KAI - not that it matters as we're working with relative index
  
  area_trans <- data.frame(ID = c(1:37),
                           transect = c("Abalario", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Hinojos",
                                        "Hinojos-Guadiamar", "Marismillas", "Muro", "Matochal", "Puntal",
                                        "RBD-este", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea",  
                                        "Sabinar-Mogea", "RBD-este", "RBD-este", "Algaida-Sotos", "Algaida-Sotos", 
                                        "Algaida-Sotos", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Coto del Rey", 
                                        "Coto del Rey", "Coto del Rey", "Coto del Rey", "Coto del Rey", "Hinojos", 
                                        "Hinojos", "Marismillas", "Marismillas", "Marismillas", "Puntal",
                                        "Puntal", "Puntal"
                           ))
  
  length_trans <- cbind(area_trans, vect_length) %>% 
    group_by(transect) %>%
    summarise(trans_length = sum(vect_length))
  
  hab_rast <- rast(hab_file)
  transects <- project(transects, crs(hab_rast))
  
  df <- c()
  
 
    for (d in c(1:nrow(dates))) {
      tryCatch({
        file = file.path(sim_output, "maps", paste0("Rabbit_Population_distribution_", dates$year[d], "_", dates$month[d], ".csv"))
        ndvi_file = file.path("data", "model_input", "maps", "ndvi_monthly_maps", 
                              paste0(dates$year[d], "_", sprintf("%02d", dates$month[d]), ".txt"))
        
        sim_mat <- as.matrix(read.csv(file, header = F))
        sim_rast <- terra::rasterize(sim_mat, hab_rast)
        values(sim_rast) <- sim_mat  
        
        ndvi_mat <- as.matrix(read.table(ndvi_file, skip = 1, header = F, sep = " "))
        ndvi_rast <- terra::rasterize(ndvi_mat, hab_rast)
        values(ndvi_rast) <- ndvi_mat
        
        #--------------------------------------------------------------------------------------------------
        # RBD nigth censuses
        #--------------------------------------------------------------------------------------------------
        date <- paste("15", dates$month[d], dates$year[d], sep = "/") 
        
        cellsT <- terra::extract(sim_rast, transects)
        cellsN <- terra::extract(ndvi_rast, transects)
        
        # calculate "detected rabbits"
        cellsD <- cellsT
        cellsD$last[which(!is.na(cellsT$last))] <- rbinom(length(which(!is.na(cellsT$last))), cellsT$last[which(!is.na(cellsT$last))],
                                                          boot::inv.logit(obs_p + (cellsN$last[which(!is.na(cellsD$last))] * obs_p_ndvi)))
        
        
        suppressMessages(
          df_yr <- left_join(cellsD, area_trans) %>%
            group_by(transect) %>%
            summarise(totR = sum(last, na.rm = T)) %>%
            ungroup() %>%
            left_join(., length_trans) %>%
            mutate(KAI = totR/trans_length) %>% 
            select(transect, KAI) %>%
            mutate(Fecha = as.Date(date, format = "%d/%m/%Y")) %>%
            filter(transect %in% c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea", "RBD-este", "Puntal", "Marismillas", "Abalario", "Hinojos"))
        )
        df <- rbind(df, df_yr)
        
      }, error = function(e) {
        # Silently catch and ignore errors
        NULL
      }, warning = function(w) {
        # Silently catch and ignore warnings
        NULL
      })
    }
  }))
  return(df)
}



csvToRaster <- function(fileName, habitat_raster, return_df = F, plot = T, 
                        save_tiff = F, tiff_name = NA) {
  
  mat_status <- as.matrix(read.csv(fileName, header = F))
  
  df <- terra::rasterize(mat_status, habitat_raster)
  values(df) <- mat_status
  
  if (return_df) {
    return(df)  
  } 
  
  if (plot) {
    plot(df)
  }
  
  if(save_tiff) {
    writeRaster(df, filename = tiff_name, overwrite = T)
  }
  
}







