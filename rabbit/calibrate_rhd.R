#
# RHD calibration, steps 0 and 1 (see backcast/backcast_report.md, "RHD disease model")
#
# Step 0: sensitivity grid over the RHD settings, to see which control the outbreak interval and
#         the declines after an arrival.
# Step 1: disease-only rejection ABC on the settings that matter, with the habitat and demography settings
#         fixed (one posterior row of the rabbit calibration).
# Step 2 (the joint recalibration) uses model_fn() with rhd_settings / rhd_param_names /
#         rhd_target in rabbit/Workflow_calibration.R.
#
# The K maps and climate files are prepared once (prepare_rhd_inputs); after that each run calls
# the Pascal program directly and takes seconds. Run from the project root.

library(terra)
library(lubridate)
library(dplyr)
library(ncdf4)

sapply(list.files(file.path("rabbit", "functions"), pattern = "\\.R$", full.names = TRUE), source)
source(file.path("backcast", "functions", "backcast_helpers.R"))

#-------------------------------------------------------------------------------
# Inputs: K maps of one posterior row and the rabbit climate files
#-------------------------------------------------------------------------------

prepare_rhd_inputs <- function(
    start_year, end_year, out_dir,
    posterior_file      = file.path("rabbit", "results", "posterior_samples.csv"),
    posterior_row       = NULL,     # NULL = the row with the highest weight
    template            = file.path("data", "environmental_data", "template_raster_500.tif"),
    temperature_files   = c(file.path("data", "environmental_data", "CDS", "kriged_temperature_v2"),
                            file.path("data", "backcast", "CDS", "kriged_temperature")),
    precipitation_files = c(file.path("data", "environmental_data", "CDS", "precipitation", "precipitation.nc"),
                            file.path("data", "backcast", "CDS", "precipitation")),
    ndvi_mean_file      = file.path("data", "environmental_data", "monthly_ndvi_v2", "mean.tif"),
    ndvi_sd_file        = file.path("data", "environmental_data", "monthly_ndvi_v2", "sd.tif"),
    flood_files         = file.path("data", "backcast", "monthly_floodmaps"),
    height_files        = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif")
){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  DT_template <- rast(template)

  window <- function(r){
    t <- as.Date(time(r))
    r[[t >= as.Date(sprintf("%d-11-01", start_year - 1)) & t <= as.Date(sprintf("%d-12-01", end_year))]]
  }
  existing <- function(paths) paths[file.exists(paths)]

  r_temp   <- window(load_monthly_stack(existing(temperature_files), DT_template, pattern = "_Kriged\\.nc$"))
  r_precip <- window(load_monthly_stack(existing(precipitation_files), DT_template, pattern = "\\.nc$"))
  start_date <- as.Date(sprintf("%d-01-01", start_year))
  r_mean <- window(fill_leading_climatology(read_monthly_ndvi(ndvi_mean_file), start_date))
  r_sd   <- window(fill_leading_climatology(read_monthly_ndvi(ndvi_sd_file), start_date))
  r_flood  <- load_and_format_rasters(flood_files, DT_template, aggregate = TRUE, fun = "max",
                                      temporal_agg = "month", temporal_fun = "max")
  r_height <- load_and_format_rasters(height_files, DT_template, aggregate = TRUE, fun = "mean")

  checked <- check_time_span_list(list(temperature = r_temp, precipitation = r_precip,
                                       meanNDVI = r_mean, sdNDVI = r_sd),
                                  start_year, end_year, fill_gaps = TRUE)

  clim_dir <- file.path(out_dir, "clim_variables")
  calculate_BM_cDM(temp = checked$temperature, K = TRUE, precip = checked$precipitation,
                   year_min = start_year, year_max = end_year, result_dir = clim_dir,
                   coord_rast = DT_template)

  post <- read.csv(posterior_file)
  p <- if (is.null(posterior_row)) post[which.max(post$weights), ] else post[posterior_row, ]
  k_file <- create_K_map_gaussian(
    Cov_list   = list(ndvi_mean = checked$meanNDVI, ndvi_sd = checked$sdNDVI, flood = r_flood, height = r_height),
    param_list = list(K = p$K_max, N_opt = p$ndvi_opt, N_sig = p$ndvi_sigma, S_opt = p$ndvi_sd_opt,
                      S_sig = p$ndvi_sd_sigma, H_opt = p$height_opt, H_sig = p$height_sigma,
                      N_w = p$ndvi_w, S_w = p$ndvi_sd_w, H_w = p$height_w),
    output_dir = out_dir, start_year = start_year, end_year = end_year, study_area = DT_template)

  inputs <- list(k_file = normalizePath(k_file), clim_dir = normalizePath(clim_dir),
                 start_year = start_year, end_year = end_year,
                 dens_opt = round(p$dens_opt), lambda = p$R_lambda, sigma = p$R_sigma)
  saveRDS(inputs, file.path(out_dir, "rhd_inputs.rds"))
  inputs
}

#-------------------------------------------------------------------------------
# One run of the Pascal program with given RHD settings
#-------------------------------------------------------------------------------

run_rhd <- function(inputs, settings, seed = -1,
                    exe   = file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"),
                    rates = file.path("data", "rabbit", "Rabbit_demography_values.txt")){
  run_dir <- tempfile("rhd_run_")
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)
  settings$RHD_log <- 1
  write_rhd_settings(settings, file.path(run_dir, "rhd_settings.txt"))

  args <- c(run_dir, inputs$start_year, inputs$end_year, normalizePath(rates),
            inputs$k_file, inputs$k_file, inputs$clim_dir, "0", "0",
            inputs$dens_opt, inputs$lambda, inputs$sigma, seed, file.path(run_dir, "rhd_settings.txt"))
  out <- system2(normalizePath(exe), args = as.character(args), stdout = TRUE, stderr = TRUE, timeout = 1800)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    stop("Model run failed: ", paste(tail(out, 5), collapse = "\n"))
  }
  read_rhd_logs(run_dir)
}

#-------------------------------------------------------------------------------
# Statistics and targets
#-------------------------------------------------------------------------------

# targets: list(interval = list(scale = "park", range = c(24, 36), from = "YYYY-MM" or NULL),
#               declines = list(list(arrival = "YYYY-MM", range = c(-0.75, -0.55)), ...))
rhd_default_targets <- function(){
  list(interval = list(scale = "park", range = c(24, 36), from = NULL),
       declines = list())
}

range_distance <- function(x, range){
  if (is.na(x)) return(5)
  width <- diff(range)
  if (x < range[1]) (range[1] - x) / width else if (x > range[2]) (x - range[2]) / width else 0
}

rhd_run_stats <- function(logs, targets){
  from <- targets$interval$from
  from <- if (is.null(from)) NULL else { a <- as.integer(strsplit(from, "-")[[1]]); rhd_month_index(a[1], a[2]) }
  iv <- rhd_intervals(logs, scale = targets$interval$scale, from = from)
  stats <- c(interval_median  = if (length(iv) > 0) median(iv) else NA,
             interval_distance = rhd_interval_distance(iv, targets$interval$range),
             cell_interval_median = median(rhd_cell_intervals(logs$outbreaks)),
             mean_rabbits = mean(logs$monthly$rabbits),
             min_rabbits  = min(logs$monthly$rabbits),
             mean_immune_share = mean(logs$monthly$immune / pmax(logs$monthly$rabbits, 1)),
             # 0 while the population persists; 5 if it dies out (no rabbits left in a month)
             persistence_distance = if (min(logs$monthly$rabbits) > 0) 0 else 5)
  for (i in seq_along(targets$declines)) {
    d <- rhd_decline(logs$monthly, targets$declines[[i]]$arrival)
    stats[paste0("decline_", i)] <- d
    stats[paste0("decline_", i, "_distance")] <- range_distance(d, targets$declines[[i]]$range)
  }
  stats
}

#-------------------------------------------------------------------------------
# Step 0: sensitivity grid
#-------------------------------------------------------------------------------

# grid: data.frame with one column per RHD setting to vary (e.g. expand.grid(...)); `base` holds
# the other settings. Each combination is run `reps` times.
rhd_sensitivity_grid <- function(inputs, grid, base = rhd_default_settings(),
                                 targets = rhd_default_targets(), reps = 3, out_file = NULL){
  res <- vector("list", nrow(grid) * reps)
  k <- 0
  for (i in seq_len(nrow(grid))) {
    settings <- modifyList(base, as.list(grid[i, , drop = FALSE]))
    for (r in seq_len(reps)) {
      k <- k + 1
      stats <- tryCatch(rhd_run_stats(run_rhd(inputs, settings, seed = 1000 * i + r), targets),
                        error = function(e) { message("Run failed: ", conditionMessage(e)); NULL })
      if (!is.null(stats)) res[[k]] <- data.frame(grid[i, , drop = FALSE], rep = r, as.list(stats))
    }
    message(sprintf("grid %d/%d done", i, nrow(grid)))
  }
  out <- bind_rows(res)
  if (!is.null(out_file)) write.csv(out, out_file, row.names = FALSE)
  out
}

#-------------------------------------------------------------------------------
# Step 1: disease-only ABC
#-------------------------------------------------------------------------------

# Rejection ABC: `n_sim` draws from the priors, each run once; the `accept_share` closest to the
# targets form the posterior. (Plain rejection rather than EasyABC: with two or three settings and
# runs of seconds it is fast enough, and EasyABC fails when a statistic has the same value in every
# run, which the persistence statistic often has.)
# priors: named list of c("unif", lower, upper); a name starting with "log10_" is sampled on the
# log10 scale (e.g. log10_RHD_spark = c("unif", -5, -2)). Distance: the root of the summed squares
# of the interval distance, the persistence distance and the decline distances (all 0 on target).
rhd_abc <- function(inputs, priors, base = rhd_default_settings(), targets = rhd_default_targets(),
                    n_sim = 300, accept_share = 0.1, out_file = NULL){
  param_names <- names(priors)
  stat_names  <- c("interval_distance", "persistence_distance",
                   if (length(targets$declines) > 0) paste0("decline_", seq_along(targets$declines), "_distance"))
  
  draws <- as.data.frame(lapply(priors, function(p){
    if (p[1] != "unif") stop("Only uniform priors are supported")
    runif(n_sim, as.numeric(p[2]), as.numeric(p[3]))
  }))
  
  results <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    settings <- base
    for (nm in param_names) {
      if (startsWith(nm, "log10_")) settings[[sub("^log10_", "", nm)]] <- 10^draws[i, nm]
      else settings[[nm]] <- draws[i, nm]
    }
    stats <- tryCatch(rhd_run_stats(run_rhd(inputs, settings, seed = i), targets),
                      error = function(e) { message("Run failed: ", conditionMessage(e)); NULL })
    results[[i]] <- if (is.null(stats)) NULL else as.data.frame(as.list(stats))
    if (i %% 25 == 0) message(sprintf("ABC: %d/%d runs", i, n_sim))
  }
  
  ok <- !vapply(results, is.null, logical(1))
  sims <- cbind(draws[ok, , drop = FALSE], bind_rows(results[ok]))
  sims$distance <- sqrt(rowSums(as.matrix(sims[, stat_names])^2))
  sims <- sims[order(sims$distance), ]
  sims$accepted <- seq_len(nrow(sims)) <= ceiling(accept_share * nrow(sims))
  
  if (!is.null(out_file)) write.csv(sims, out_file, row.names = FALSE)
  list(posterior = sims[sims$accepted, ], all = sims)
}

#-------------------------------------------------------------------------------
# Example (not run when sourced)
#-------------------------------------------------------------------------------

if (FALSE) {
  # Endemic period with census data; RHD present from the start
  inputs <- prepare_rhd_inputs(2005, 2024, file.path("rabbit", "results", "rhd_calibration", "inputs"))

  grid <- expand.grid(RHD_f_crit = c(0.3, 0.5, 0.7), RHD_spark = c(0.0005, 0.002, 0.01),
                      RHD_spread = c(0, 0.2, 0.5), RHD_mortality = c(0.4, 0.7))
  sens <- rhd_sensitivity_grid(inputs, grid, reps = 3,
                               out_file = file.path("rabbit", "results", "rhd_calibration", "sensitivity_grid.csv"))

  abc <- rhd_abc(inputs,
                 priors = list(RHD_f_crit = c("unif", 0.2, 0.9),
                               log10_RHD_spark = c("unif", -4, -1.5),
                               RHD_spread = c("unif", 0, 1)),
                 n_sim = 300, accept_share = 0.1,
                 out_file = file.path("rabbit", "results", "rhd_calibration", "rhd_abc_runs.csv"))

  # Phase 2 period: both arrivals, with the two documented declines as extra targets
  # targets <- modifyList(rhd_default_targets(),
  #   list(declines = list(list(arrival = "1989-06", range = c(-0.75, -0.55)),
  #                        list(arrival = "2011-11", range = c(-1, -0.8)))))
}
