
# Try to recover a timestep from a file name when the raster/file itself
# doesn't carry time information. Handles daily (e.g. 2023-01-15,
# 2023_01_15, 20230115), monthly (2023-01, 2023_01) and yearly (2023)
# naming conventions, checked in that order of specificity.
# Kept identical to the copies in NDVI_model_component/code/create_model_dataframe.R
# and NDVI_forecast_component/code/create_model_dataframe.R - update all three together.
parse_time_from_filename <- function(f){
  bn <- tools::file_path_sans_ext(basename(f))

  # daily
  m <- regmatches(bn, regexpr("(\\d{4})[-_]?(\\d{2})[-_]?(\\d{2})(?!\\d)", bn, perl = TRUE))
  if(length(m) == 1 && nzchar(m)){
    parts <- regmatches(m, regexec("(\\d{4})[-_]?(\\d{2})[-_]?(\\d{2})", m))[[1]]
    d <- suppressWarnings(as.Date(sprintf("%s-%s-%s", parts[2], parts[3], parts[4])))
    if(!is.na(d)) return(d)
  }

  # monthly
  m <- regmatches(bn, regexpr("(\\d{4})[-_](\\d{2})(?!\\d)", bn, perl = TRUE))
  if(length(m) == 1 && nzchar(m)){
    parts <- regmatches(m, regexec("(\\d{4})[-_](\\d{2})", m))[[1]]
    d <- suppressWarnings(as.Date(sprintf("%s-%s-01", parts[2], parts[3])))
    if(!is.na(d)) return(d)
  }

  # yearly
  m <- regmatches(bn, regexpr("(?<!\\d)(19|20)\\d{2}(?!\\d)", bn, perl = TRUE))
  if(length(m) == 1 && nzchar(m)){
    d <- suppressWarnings(as.Date(sprintf("%s-01-01", m)))
    if(!is.na(d)) return(d)
  }

  NULL
}

load_and_format_rasters <- function(files, template_raster,
                                    aggregate = F, fun = "sd",
                                    temporal_agg = NULL, temporal_fun = "mean"){
  if(length(files) == 1){
    if(file_test("-d", files)){
      files = list.files(files, full.names = T)
    }
  }

  if(length(files) == 0){
    stop('no files found')
  }

  r_list <- list()

  for(f in files){
    if(tools::file_ext(f) %in% c("nc")){

      temp_r <- rast(f)

      if(is.na(crs(temp_r)) || crs(temp_r) == ""){
        crs(temp_r) <- crs(template_raster)
      }

      nc <- nc_open(f)

      # The time dimension isn't always called "time" (e.g. CDS/ERA5 files)
      dim_names <- names(nc$dim)
      time_dim  <- dim_names[grepl("^time$", dim_names, ignore.case = TRUE)]
      if(length(time_dim) == 0){
        time_dim <- dim_names[grepl("time", dim_names, ignore.case = TRUE)]
      }

      if(length(time_dim) > 0){
        time_dim <- time_dim[1]

        time_values <- ncvar_get(nc, time_dim)
        time_units  <- ncatt_get(nc, time_dim, "units")$value

        nc_close(nc)

        # Extract origin dynamically
        origin <- sub("^.* since", "", time_units)

        time_unit_type <- regmatches(time_units, regexpr("^.* since", time_units))
        time_unit_type <- sub(" since", "", time_unit_type)

        if(!(time_unit_type %in% c("milliseconds", "seconds"))){
          stop('Time unit of the raster is not milli(seconds) since ....')
        }

        if(time_unit_type == "milliseconds"){
          time_values <- time_values / 1000
        }

        # Convert milliseconds → seconds and decode
        times <- as.POSIXct(
          time_values,
          origin = origin,
          tz = "UTC"
        )
        time(temp_r) <- times
      } else {
        nc_close(nc)
      }

    } else {
      temp_r <- rast(f)
    }

    # Some files (no "time" dimension in the nc, or non-nc formats like
    # GeoTIFF that never store time) don't carry time information at all,
    # but encode the timestep - daily, monthly or yearly - in the file name.
    if(all(is.na(time(temp_r)))){
      file_time <- parse_time_from_filename(f)
      if(is.null(file_time)){
        message("no time information found (no time dimension/metadata, and could not parse a date from the file name) in ", f)
      }
      time(temp_r) <- rep(file_time, nlyr(temp_r))
    }

    r_list[[length(r_list) + 1]] <- temp_r

  }

  r <- if(length(r_list) == 1) r_list[[1]] else rast(r_list)

  if(!same.crs(r, template_raster)){
    r <- project(r, crs(template_raster))
  }

  if(!is.null(temporal_agg)){
    temporal_agg <- match.arg(temporal_agg, c("month", "year"))
    r_time <- time(r)
    if(all(is.na(r_time))){
      stop("Cannot temporally aggregate: no time information available on the rasters.")
    }
    period <- floor_date(as.Date(r_time), temporal_agg)
    r <- terra::tapp(r, index = period, fun = temporal_fun, na.rm = TRUE)
    time(r) <- sort(unique(period))
  }


  if(aggregate){
    res_r = res(r)
    res_temp = res(template_raster)
    r <- terra::aggregate(r, fact = res_temp/res_r, fun = fun, na.rm = T)
  }
  
  r <- resample(r, template_raster)

  gc()
  suppressWarnings(terra::tmpFiles(remove = TRUE))

  return(r)
}


extract_year_month <- function(path) {
  bn <- basename(path)
  yr <- as.integer(sub("^(\\d{4})_(\\d{1,2}).*$", "\\1", bn))
  mo <- as.integer(sub("^(\\d{4})_(\\d{1,2}).*$", "\\2", bn))
  c(yr, mo)
}


parse_time_from_names <- function(rast_obj) {
  names_vec <- names(rast_obj)
  parsed <- suppressWarnings(as.Date(gsub(".*?(\\d{4})_(\\d{1,2}).*", "\\1-\\2-01", names_vec)))
  if (all(is.na(parsed))) return(NULL)
  parsed
}

get_stack_time <- function(rast_obj, time_file = NULL) {
  if (!is.null(time_file) && file.exists(time_file)) {
    time_vals <- readRDS(time_file)
    if (is.numeric(time_vals) && all(time_vals %in% 1:12)) {
      return(as.integer(time_vals))
    }
    if (inherits(time_vals, "Date")) {
      return(as.POSIXct(time_vals))
    }
    if (inherits(time_vals, "POSIXct")) {
      return(time_vals)
    }
    return(as.POSIXct(time_vals, origin = "1970-01-01", tz = "UTC"))
  }
  time_vals <- time(rast_obj)
  if (all(is.na(time_vals))) {
    parsed <- parse_time_from_names(rast_obj)
    if (!is.null(parsed)) {
      return(as.POSIXct(parsed))
    }
    parsed_months <- suppressWarnings(as.integer(gsub(".*m_(\\d{1,2}).*", "\\1", names(rast_obj))))
    if (!all(is.na(parsed_months))) {
      return(as.integer(parsed_months))
    }
  }
  if (all(is.na(time_vals))) return(NULL)
  time_vals
}

parse_month_indices <- function(rast_obj, time_file = NULL) {
  if (!is.null(time_file) && file.exists(time_file)) {
    time_vals <- readRDS(time_file)
    if (is.numeric(time_vals) && all(time_vals %in% 1:12)) {
      return(as.integer(time_vals))
    }
  }
  parsed <- suppressWarnings(as.integer(gsub(".*m_(\\d{1,2}).*", "\\1", names(rast_obj))))
  if (all(is.na(parsed))) stop("Unable to determine scenario month indices")
  parsed
}

build_variant_climate_stacks <- function(temp_obs_files,
                                         precip_obs_file,
                                         temp_scen_file,
                                         temp_scen_time,
                                         precip_scen_file,
                                         precip_scen_time,
                                         variant,
                                         scenario,
                                         output_dir,
                                         end_year) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  temp_obs <- rast(temp_obs_files)
  precip_obs <- rast(precip_obs_file)
  
  temp_obs <- project(temp_obs, crs(DT_template))
  precip_obs <- project(precip_obs, crs(DT_template))
  
  temp_obs <- resample(temp_obs, DT_template) - 273.15
  precip_obs <- resample(precip_obs, DT_template)
  
  temp_obs_time <- get_stack_time(temp_obs)
  precip_obs_time <- get_stack_time(precip_obs)
  
  if (is.null(temp_obs_time) || is.null(precip_obs_time)) {
    stop("Unable to determine observed climate time information")
  }
  time(temp_obs) <- temp_obs_time
  time(precip_obs) <- precip_obs_time
  
  observed_dates <- sort(intersect(floor_date(as.Date(time(temp_obs)), "month"), 
                                   floor_date(as.Date(time(precip_obs)), "month")))
  if (length(observed_dates) == 0) stop("No common observed climate dates available")
  
  end_date <- as.Date(paste0(end_year, "-12-01"))
  if (variant == "observed_until_2024") {
    observed_keep <- observed_dates[observed_dates <= as.Date("2024-12-01")]
  } else {
    observed_keep <- observed_dates[observed_dates <= end_date]
  }
  if (length(observed_keep) == 0) {
    stop("No observed data selected for variant ", variant)
  }
  
  temp_obs_stack <- temp_obs[[floor_date(as.Date(time(temp_obs)), "month") %in% observed_keep]]
  precip_obs_stack <- precip_obs[[floor_date(as.Date(time(precip_obs)), "month") %in% observed_keep]]
  
  max_obs_date <- max(observed_keep)
  target_dates <- seq(max_obs_date + months(1), end_date, by = "1 month")
  if (length(target_dates) > 0) {
    temp_scen <- rast(temp_scen_file)
    precip_scen <- rast(precip_scen_file)
    temp_scen_months <- parse_month_indices(temp_scen, temp_scen_time)
    precip_scen_months <- parse_month_indices(precip_scen, precip_scen_time)
    if (length(temp_scen_months) != nlyr(temp_scen)) stop("Temperature scenario raster layer count and month metadata mismatch")
    if (length(precip_scen_months) != nlyr(precip_scen)) stop("Precipitation scenario raster layer count and month metadata mismatch")
    
    temp_future <- rast(lapply(target_dates, function(dt) {
      m <- month(dt)
      idx <- which(temp_scen_months == m)
      if (length(idx) != 1) stop("Unable to find scenario temperature layer for month ", m)
      subset(temp_scen, idx)
    }))
    precip_future <- rast(lapply(target_dates, function(dt) {
      m <- month(dt)
      idx <- which(precip_scen_months == m)
      if (length(idx) != 1) stop("Unable to find scenario precipitation layer for month ", m)
      subset(precip_scen, idx)
    }))
    time(temp_future) <- as.POSIXct(target_dates)
    time(precip_future) <- as.POSIXct(target_dates)
    
    temp_comb <- c(temp_obs_stack, temp_future)
    precip_comb <- c(precip_obs_stack, precip_future)
    time(temp_comb) <- c(time(temp_obs_stack), time(temp_future))
    time(precip_comb) <- c(time(precip_obs_stack), time(precip_future))
  } else {
    temp_comb <- temp_obs_stack
    precip_comb <- precip_obs_stack
    time(temp_comb) <- time(temp_obs_stack)
    time(precip_comb) <- time(precip_obs_stack)
  }
  
  temp_file <- file.path(output_dir, paste0("temp_", variant, "_", scenario, ".tif"))
  precip_file <- file.path(output_dir, paste0("precip_", variant, "_", scenario, ".tif"))
  writeRaster(temp_comb, temp_file, overwrite = TRUE)
  writeRaster(precip_comb, precip_file, overwrite = TRUE)
  
  list(temp = temp_file, precip = precip_file)
}

save_population_mean_rasters <- function(pop_maps, out_dir, variant, scenario) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(pop_maps)) {
    mean_rast <- pop_maps[[name]]$mean
    if (is.null(mean_rast)) next
    out_file <- file.path(out_dir, paste0(variant, "_", scenario, "_", name, "_mean.tif"))
    writeRaster(mean_rast, out_file, overwrite = TRUE)
  }
}



# Check if the time span of a raster covers the full year range given 
# + November and December before the start year, for the rabbit variables
# does not check timesteps (months, days etc.)

check_time_span <- function(r, start_year, end_year, raster_name = NULL) {
  t <- time(r)
  
  label <- if (!is.null(raster_name)) raster_name else deparse(substitute(r))
  
  if (is.null(t) || all(is.na(t))) {
    stop(sprintf("'%s' has no time information set (time(r) is NULL/NA).", label))
  }
  
  expected_start <- as.Date(paste0(start_year, "-01-01"))
  expected_end   <- as.Date(paste0(end_year, "-10-31"))
  
  actual_start <- min(t, na.rm = TRUE)
  actual_end   <- max(t, na.rm = TRUE)
  
  if (inherits(actual_start, "POSIXct")) {
    actual_start <- as.Date(actual_start)
    actual_end   <- as.Date(actual_end)
  } else if (is.numeric(actual_start)) {
    actual_start <- as.Date(paste0(actual_start, "-01-01"))
    actual_end   <- as.Date(paste0(actual_end, "-12-31"))
  }
  
  ok_start <- actual_start <= expected_start
  ok_end   <- actual_end >= expected_end
  
  if (!ok_start || !ok_end) {
    msg <- sprintf(
      "'%s' time span (%s to %s) does not cover required range (%s to %s).",
      label, actual_start, actual_end, expected_start, expected_end
    )
    if (!ok_start) msg <- paste(msg, "Missing coverage at START.")
    if (!ok_end)   msg <- paste(msg, "Missing coverage at END.")
    stop(msg, call. = FALSE)
  }
  
  invisible(TRUE)
}

check_time_span_list <- function(raster_list, start_year, end_year) {
  names_list <- names(raster_list)
  if (is.null(names_list)) names_list <- paste0("raster_", seq_along(raster_list))
  
  for (i in seq_along(raster_list)) {
    check_time_span(raster_list[[i]], start_year, end_year, raster_name = names_list[i])
  }
  invisible(TRUE)
}
