
create_model_rasters <- function(template, ndvi_file, krig_dir, precip_file){
  
  # NDVI data --------------------------------------------------------------------
  cat('using uncleaned ndvi data right now!!!')
  
  ndvi_stack <- rast(ndvi_file)
  ndvi_stack <- project(ndvi_stack, crs(template))
  ndvi_stack <- resample(ndvi_stack, template, method = "bilinear")
  
  ym <- format(as.Date(names(ndvi_stack)), "%Y-%m")
  ndvi_monthly_max <- tapp(ndvi_stack, index = ym, fun = max, na.rm = TRUE)
  date <- gsub("X", "", names(ndvi_monthly_max))
  date <- as.Date(paste0(gsub("\\.", "-", date), "-15"), format = "%Y-%m-%d")
  terra::time(ndvi_monthly_max) <- date
  
  # NDVI neighbours (for spatial autocorrelation) --------------------------------
  
  ndvi_auto_max3 <- focal(
    ndvi_stack,
    w =  3,
    fun = max,
    na.rm = TRUE
  )
  terra::time(ndvi_auto_max3) <- as.Date(names(ndvi_auto_max3))
  
  # temperature-------------------------------------------------------------------
  
  temp_files = list.files(krig_dir, pattern = "Kriged.nc$", full.names = T)
  temp_stack = lapply(temp_files, rast) %>% rast(.)
  temp_stack = project(temp_stack, crs(template))
  temp_stack = resample(temp_stack, template, method = "bilinear")
  temp_stack = temp_stack - 273.15
  # names(temp_stack) = rep("temp", length(names(temp_stack)))
  
  # Precipitation -------------------------------------------------------------------
  
  precip_stack = rast(precip_file)
  precip_stack = project(precip_stack, crs(template))
  precip_stack = resample(precip_stack, template, method = "bilinear")
  
  # transform precipitation from mean m/day to total mm/month
  dates <- time(precip_stack)
  scaling <- 1000 * days_in_month(dates)
  precip_scaled <- precip_stack
  
  for (i in 1:nlyr(precip_stack)) {
    precip_scaled[[i]] <- precip_stack[[i]] * scaling[i]
  }
  
  precip_scaled <- focal(precip_scaled, w=15, fun=mean, NAonly=T, na.rm=T)
  # names(precip_scaled) <- rep("precip", length(names(precip_scaled)))
  
  
  return(list("ndvi_stack" = ndvi_monthly_max,
              "ndvi_neighbour" = ndvi_auto_max3,
              "temp_stack" = temp_stack,
              "precip_stack" = precip_scaled))
}


create_model_dataframe <- function(ndvi_stack, ndvi_auto_max3,
                                   temp_stack, precip_stack) {
  
  #-------------------------------------------------------------------------------
  # Format data  
  #-------------------------------------------------------------------------------
  
  ndvi_data <- as.data.frame(ndvi_stack, xy = TRUE, time = TRUE, wide = FALSE) %>%
    rename(ndvi = values)
  ndvi_data$layer = NULL
  ndvi_data$year <- year(ndvi_data$time)
  ndvi_data$month <- month(ndvi_data$time)
  ndvi_data <- ndvi_data %>% 
    group_by(x, y, year, month) %>% 
    summarise(max_ndvi = max(ndvi, na.rm = T))
  
  # spatial autocorrelation
  ndvi_auto_data <- as.data.frame(ndvi_auto_max3, xy = TRUE, time = TRUE, wide = FALSE) %>%
    rename(ndvi_max3 = values)
  ndvi_auto_data$year <- year(ndvi_auto_data$time)
  ndvi_auto_data$month <- month(ndvi_auto_data$time)
  ndvi_auto_data$layer <- NULL
  ndvi_auto_data$time <- NULL
  
  
  # Temperature
  temp_data <- as.data.frame(temp_stack, xy = TRUE, time = TRUE, wide = FALSE) %>%
    rename(temp = values)
  temp_data$year <- year(temp_data$time)
  temp_data$month <- month(temp_data$time)
  temp_data$layer <- NULL
  temp_data$time <- NULL
  
  
  # Precipitation
  precip_data <- as.data.frame(precip_stack, xy = TRUE, time = TRUE, wide = FALSE) %>%
    rename(precip = values)
  precip_data$year <- year(precip_data$time)
  precip_data$month <- month(precip_data$time)
  precip_data$layer <- NULL
  precip_data$time <- NULL
  
  # Create main dataframe --------------------------------------------------------
  data_all <- left_join(ndvi_data, ndvi_auto_data, by = c('x', 'y', 'year', 'month'))
  data_all <- left_join(data_all, temp_data, by = c('x', 'y', 'year', 'month'))
  data_all <- left_join(data_all, precip_data, by = c('x', 'y', 'year', 'month'))
  data_all <- na.omit(data_all)
  data_all <- data_all[with(data_all, order(year, month)), ]
  
  # create pixel ID for the next steps
  setDT(data_all)
  data_all[, pixel_id := .GRP, by = .(x, y)]
  
  # temporal autocorrelation
  data_all[, ndvi_lag1 := shift(max_ndvi, 1), by = pixel_id]
  data_all[, ndvi_lag3 := shift(max_ndvi, 3), by = pixel_id]
  data_all[, ndvi_lag12 := shift(max_ndvi, 12), by = pixel_id]
  
  # spatial auto correlation
  data_all[, ndvi_max3_lag1 := shift(ndvi_max3, 1), by = pixel_id]
  
  # seasonal encoding
  data_all[, month_sin := sin(2 * pi * month / 12)]
  data_all[, month_cos := cos(2 * pi * month / 12)]
  
  # Keep only complete rows
  data_all <- na.omit(data_all)
  
  return(data_all)
}