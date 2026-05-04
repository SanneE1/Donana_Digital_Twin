library(terra)
library(lubridate)

load_stack_with_time <- function(tif_path, time_path) {
  r <- rast(tif_path)
  t <- readRDS(time_path)
  stopifnot(nlyr(r) == length(t))
  terra::time(r) <- t
  return(r)
}


forecast_ndvi_raster_multivar <- function(
    model_dir,
    climate_type = "years",
    climate_years = c(2015, 2016),
    horizon = 24,
    stochastic = FALSE,
    residuals = NULL
) {
  
  if (!climate_type %in% c("years", "mean")) {
    stop("climate_type can only be 'years' or 'mean'")
  }
  
  if (climate_type == "years") {
    if (horizon != length(climate_years) * 12) {
      stop("There's a mismatch between the number of months to forecast (horizon) and the number of climate_years provided")
    }
  }
  
  model <- xgb.load(file.path(model_dir, "ndvi_xgb_model.json"))
  
  ndvi_stack <- load_stack_with_time(tif_path = file.path(model_dir, "ndvi_stack.tif"),
                                     time_path = file.path(model_dir, "ndvi_time.rds"))
  temp_stack <- load_stack_with_time(tif_path = file.path(model_dir, "temp_stack.tif"),
                                     time_path = file.path(model_dir, "temp_time.rds"))
  precip_stack <- load_stack_with_time(tif_path = file.path(model_dir, "precip_stack.tif"),
                                     time_path = file.path(model_dir, "precip_time.rds"))
  
  climate_stacks <- list(
    precip = precip_stack,
    temp = temp_stack
  )
  
  features <- c(
    "ndvi_lag1", "ndvi_lag3", "ndvi_lag12",
    "month_sin", "month_cos",
    "ndvi_max3_lag1",
    "temp", "precip"
  )
  
  if (stochastic) {
    residuals <- readRDS(file.path(output_dir, "residuals.rds"))
  }
  
  # --- checks ---
  stopifnot(is.list(climate_stacks))
  stopifnot(all(names(climate_stacks) %in% features))
  
  # time handling
  dates <- time(ndvi_stack)
  last_date <- max(dates)
  
  results <- vector("list", 12 + horizon)
  
  for(i in 1:12) { 
    results[[i]] <- ndvi_stack[[nlyr(ndvi_stack) - (12-i)]]  
  }
  
  
  for (l in 13:length(results)) {
    
    h = l - 12
    
    new_date <- last_date %m+% months(h)
    m <- month(new_date)
    
    # initialize lag rasters
    ndvi_lag1  <- ndvi_stack[[l-1]]
    ndvi_lag3  <- ndvi_stack[[l-3]]
    ndvi_lag12 <- ndvi_stack[[l-12]]
    
    # spatial lag from LAST observed NDVI only
    ndvi_max3_lag1 <- focal(
      ndvi_lag1,
      w = 3,
      fun = max,
      na.rm = TRUE
    )
    
    # --- seasonal features ---
    month_sin <- ndvi_lag1
    values(month_sin) <- sin(2 * pi * m / 12)
    
    month_cos <- ndvi_lag1
    values(month_cos) <- cos(2 * pi * m / 12)
    
    # --- climate extraction ---
    ref_year <- climate_years[(h - 1) %% length(climate_years) + 1]
    
    clim_layers <- list()
    
    for (var in names(climate_stacks)) {
      
      stack <- climate_stacks[[var]]
      tvec  <- time(stack)
      
      if (climate_type == "years") {
        
        ref_year <- climate_years[(h - 1) %% length(climate_years) + 1]
        idx <- which(year(tvec) == ref_year & month(tvec) == m)
        
        if (length(idx) == 0) {
          stop(paste("No climate data for", var, ref_year, m))
        }
        
        layer <- stack[[idx]]
        
      } else if (climate_type == "mean") {
        
        idx <- which(month(tvec) == m)
        
        if (length(idx) == 0) {
          stop(paste("No climate data for", var, "month", m))
        }
        
        layer <- mean(stack[[idx]])
      }
      
      names(layer) <- var
      clim_layers[[var]] <- layer
    }
    
    
    clim_stack <- rast(clim_layers)
    
    # --- build feature stack explicitly ---
    X <- c(
      ndvi_lag1,
      ndvi_lag3,
      ndvi_lag12,
      ndvi_max3_lag1,
      month_sin,
      month_cos,
      clim_stack
    )
    
    names(X) <- c(
      "ndvi_lag1",
      "ndvi_lag3",
      "ndvi_lag12",
      "ndvi_max3_lag1",
      "month_sin",
      "month_cos",
      names(climate_stacks)
    )
    
    # ensure correct feature order
    X <- X[[features]]
    
    # --- prediction ---
    pred <- predict(
      X,
      model,
      fun = function(m, d) predict(m, as.matrix(d))
    )
    
    # --- stochastic option ---
    if (stochastic) {
      noise <- sample(residuals, ncell(pred), replace = TRUE)
      pred <- pred + setValues(pred, noise)
    }
    
    results[[l]] <- pred
    
    
  }
  
  # --- assemble output ---
  future_stack <- rast(results)
  
  terra::time(future_stack) <- seq(
    last_date %m-% months(12),
    by = "month",
    length.out = length(results)
  )
  
  return(future_stack)
}