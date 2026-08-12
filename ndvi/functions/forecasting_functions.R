library(terra)
library(lubridate)

load_stack_with_time <- function(tif_path, time_path) {
  r <- rast(tif_path)
  t <- readRDS(time_path)
  stopifnot(nlyr(r) == length(t))
  terra::time(r) <- t
  return(r)
}


forecast_ndvi_raster <- function(
    data_dir,
    model_dir,
    future_climate,   # named list: list(temp = SpatRaster, precip = SpatRaster)
    horizon = 24,
    stochastic = FALSE,
    residuals = NULL
) {
  
  model <- xgb.load(here::here(file.path(model_dir, "ndvi_xgb_model.json")))
  
  ndvi_stack <- load_stack_with_time(tif_path = here::here(file.path(data_dir, "ndvi_stack.tif")),
                                     time_path = here::here(file.path(data_dir, "ndvi_time.rds")))
  
  features <- c(
    "ndvi_lag1", "ndvi_lag3", "ndvi_lag12",
    "month_sin", "month_cos",
    "ndvi_max3_lag1",
    "temp", "precip"
  )
  
  if (stochastic) {
    residuals <- readRDS(here::here(file.path(model_dir, "residuals.rds")))
  }
  
  
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
    # Pull the h-th layer from the pre-built future climate rasters
    temp_layer   <- future_climate$temp[[h]]
    precip_layer <- future_climate$precip[[h]]
    names(temp_layer)   <- "temp"
    names(precip_layer) <- "precip"
    
    # --- build feature stack explicitly ---
    X <- c(
      ndvi_lag1,
      ndvi_lag3,
      ndvi_lag12,
      ndvi_max3_lag1,
      month_sin,
      month_cos,
      temp_layer,
      precip_layer
    )
    
    names(X) <- features
    
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


update_baseline <- function(future_climate, ndvi_last_date,
                            data_dir, model_dir, pred_dir,
                            baseline_pred, baseline_time
) {
  
  temp_map   <- c("warm" = "warm", "average" = "mean", "cold" = "cold")
  precip_map <- c("wet" = "wet",  "average" = "mean", "dry" = "dry")
  
  future_months <- seq(ndvi_last_date %m+% months(1), by = "month", length.out = 24)
  
  temp_layers <- lapply(seq_len(24), function(h) {
    month_num <- month(future_months[h])
    stack     <- future_climate$temp[["mean"]]
    stack[[which(time(stack) == month_num)]]
  })
  
  precip_layers <- lapply(seq_len(24), function(h) {
    month_num <- month(future_months[h])
    stack     <- future_climate$precip[["mean"]]
    stack[[which(time(stack) == month_num)]]
  })
  
  baseline_climate <- list(
    temp   = rast(temp_layers),
    precip = rast(precip_layers)
  )
  
  pred_reps <- lapply(as.list(1:5), function(x) {
    forecast_ndvi_raster(
      data_dir = data_dir,
      model_dir = model_dir, 
      future_climate = baseline_climate,
      stochastic = TRUE
    )
  }) %>% rast(.) 
  
  pred_reps <- tapp(pred_reps, "days", "mean")
  
  # Save prediction if this was a clean preset 
  stem      <- paste0("NDVI_preset_all_avg")
  tif_path  <- file.path(pred_dir, paste0(stem, ".tif"))
  time_path <- file.path(pred_dir, paste0(stem, "_time.rds"))
  
  times_to_save <- terra::time(pred_reps)
  attr(times_to_save, "ndvi_last_date") <- ndvi_last_date
  
  writeRaster(pred_reps, tif_path, overwrite = TRUE)
  saveRDS(times_to_save, time_path)
}




shiny_forecast_function <- function(temp_choices,    # character vector length 24: "warm"/"average"/"cold"
                                    precip_choices,  # character vector length 24: "wet"/"average"/"dry"
                                    future_climate,
                                    ndvi_last_date,
                                    data_dir,
                                    model_dir = file.path(data_dir, "model_info"),
                                    pred_dir = file.path(data_dir, "predictions"),
                                    baseline_pred = here::here(pred_dir, "NDVI_preset_all_avg.tif"),
                                    baseline_time = here::here(pred_dir, "NDVI_preset_all_avg_time.rds")
){
  
  if(!dir.exists(pred_dir)) {dir.create(pred_dir)}
  
  # Check timeline of saved baseline predictions and rerun if missing or up-to-date
  rerun_baseline = F
  if(file.exists(baseline_time)) {
    baseline_times <- readRDS(baseline_time)
    baseline_last  <- attr(baseline_times, "ndvi_last_date")
    if(baseline_last < ndvi_last_date) { rerun_baseline <- TRUE } 
  }
  if(!file.exists(baseline_pred) | !file.exists(baseline_time)) { rerun_baseline <- TRUE }
  
  if (rerun_baseline) { update_baseline(future_climate = future_climate, 
                                        ndvi_last_date = ndvi_last_date, 
                                        data_dir = data_dir,
                                        model_dir = model_dir, 
                                        pred_dir = pred_dir, 
                                        baseline_pred =  baseline_pred, 
                                        baseline_time = baseline_time
  )}
  
  pred_reps <- NULL
  
  # See if any of the standard presets were chosen
  preset_id <- NULL
  if (length(unique(temp_choices)) == 1 && length(unique(precip_choices)) == 1) {
    if(all(temp_choices == "average")) {
      if (all(precip_choices == "average")) {preset_id = "preset_all_avg"}
      if (all(precip_choices == "dry")) {preset_id = "preset_all_dry"}
      if (all(precip_choices == "wet")) {preset_id = "preset_all_wet"}
    }
    if(all(precip_choices == "average")) {
      if (all(temp_choices == "average")) {preset_id = "preset_all_avg"}
      if (all(temp_choices == "warm")) {preset_id = "preset_all_warm"}
      if (all(temp_choices == "cold")) {preset_id = "preset_all_cold"}
    }
    if(all(precip_choices == "dry") && all(temp_choices == "warm")) {preset_id = "preset_warm_dry"}
    if(all(precip_choices == "wet") && all(temp_choices == "cold")) {preset_id = "preset_cold_wet"}
  }
  
  # Check if preset option has been chosen, and if the saved file has the right date range
  if (!is.null(preset_id)) {
    stem      <- paste0("NDVI_", preset_id)
    tif_path  <- file.path(pred_dir, paste0(stem, ".tif"))
    time_path <- file.path(pred_dir, paste0(stem, "_time.rds"))
    
    if (file.exists(tif_path) && file.exists(time_path)) {
      cached_times <- readRDS(time_path)
      cached_last  <- attr(cached_times, "ndvi_last_date")
      
      if (!is.null(cached_last) && !is.null(ndvi_last_date) && cached_last == ndvi_last_date) {
        cat("Loading cached forecast for preset:", preset_id, "\n")
        pred_reps <- load_stack_with_time(tif_path = tif_path, time_path = time_path)
        
      } else {
        cat("Cache found but NDVI data has been updated — re-running simulation.\n")
      }
    }
  }
  
  
  if (is.null(pred_reps)) {
    temp_map   <- c("warm" = "warm", "average" = "mean", "cold" = "cold")
    precip_map <- c("wet" = "wet",  "average" = "mean", "dry" = "dry")
    
    future_months <- seq(ndvi_last_date %m+% months(1), by = "month", length.out = 24)
    
    temp_layers <- lapply(seq_len(24), function(h) {
      choice    <- temp_map[temp_choices[h]]
      month_num <- month(future_months[h])
      stack     <- future_climate$temp[[choice]]
      # layers are named m_1, m_2, etc.
      stack[[which(time(stack) == month_num)]]
    })
    
    precip_layers <- lapply(seq_len(24), function(h) {
      choice    <- precip_map[precip_choices[h]]
      month_num <- month(future_months[h])
      stack     <- future_climate$precip[[choice]]
      stack[[which(time(stack) == month_num)]]
    })
    
    future_climate <- list(
      temp   = rast(temp_layers),
      precip = rast(precip_layers)
    )
    
    pred_reps <- lapply(as.list(1:5), function(x) {
      forecast_ndvi_raster(
        data_dir = data_dir,
        model_dir = model_dir, 
        future_climate = future_climate,
        stochastic = TRUE
      )
    }) %>% rast(.) 
    
    pred_reps <- tapp(pred_reps, "days", "mean")
    
    # Save prediction if this was a clean preset 
    if (!is.null(preset_id)) {
      stem      <- paste0("NDVI_", preset_id)
      tif_path  <- file.path(pred_dir, paste0(stem, ".tif"))
      time_path <- file.path(pred_dir, paste0(stem, "_time.rds"))
      
      times_to_save <- terra::time(pred_reps)
      attr(times_to_save, "ndvi_last_date") <- ndvi_last_date
      
      writeRaster(pred_reps, tif_path, overwrite = TRUE)
      saveRDS(times_to_save, time_path)
      
      cat("Forecast cached for preset:", preset_id, "\n")
    }
  }
  
  # Plot
  ndvi_forecast <- load_stack_with_time(tif_path = baseline_pred, baseline_time)
  ndvi_mean <- global(ndvi_forecast, fun = "mean", na.rm = TRUE)
  ndvi_mean$time <- time(ndvi_forecast)
  ndvi_mean$type <- "mean"
  
  
  ndvi_current <- global(pred_reps, fun = "mean", na.rm = TRUE)
  ndvi_current$time <- time(pred_reps)
  ndvi_current$type <- "User Input"
  
  pred_df <- bind_rows(ndvi_mean, ndvi_current) 
  pred_df$type <- factor(pred_df$type, levels = c("User Input", "mean"))
  
  mean_plot <- ggplot() +
    annotate("rect",
             xmin = min(pred_df$time),
             xmax = min(pred_df$time) + months(12),
             ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 1) +
    annotate("text",
             x = min(pred_df$time) + months(6),
             y = Inf,
             label = "Observed",
             vjust = 1.5, size = 4, colour = "grey40") +
    annotate("text",
             x = min(pred_df$time) + months(24),
             y = Inf,
             label = "Predicted",
             vjust = 1.5, size = 4, colour = "grey40") +
    geom_line(data = pred_df, aes(x = time, y = mean, group = type, colour = type), linewidth = 1) +
    scale_colour_manual(name = "Climate",
                        breaks = c("mean", "User Input"),
                        values = c("black", "red")) +
    scale_x_date(expand = c(0, 0)) +
    ylab("Mean NDVI") + xlab("Date") +
    theme_minimal() + theme(text = element_text(size = 14),
                            axis.title = element_blank(),
                            plot.margin = unit(c(0,0,0,0), "pt"))
  
 
  return(list("mean_plot" = mean_plot,
              "mean_raster" = pred_reps))
  
}
