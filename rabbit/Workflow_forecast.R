
# Rabbit workflow forecasting

function(
    template                 = file.path("data", "environmental_data", "template_raster_500.tif"),
    years_range              = NULL,
    obs_transects            = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml"),
    obs_df                   = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"),
    # spatial_model            = file.path("data", "rabbit", "")
    rabbit_demographic_rates = file.path("rabbit", "output_data", "Rabbit_demography_values.txt"),
    demograph_mod            = file.path("rabbit", "output_data", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"),
    
    temperature_files        = file.path("data", "environmental_data", "CDS", "kriged_temperature"),
    is_Kelvin                = T,
    precipitation_files      = file.path("data", "environmental_data", "CDS", "precipitation.nc"),
    ndvi_files               = list.files(file.path("data", "environmental_data", "NDVI"), full.names = T),
    flood_files              = file.path("data", "environmental_data", "LAST_monthly_floodmaps"),
    height_files             = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif"),
    
    out_dir                  = file.path("rabbit", "results"),
    N_sim_samples            = 50
    
){
  
  # template                 = file.path("data", "environmental_data", "template_raster_500.tif")
  # years_range              = c(2023, 2024)
  # obs_transects            = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml")
  # obs_df                   = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv")
  # rabbit_demographic_rates = file.path("rabbit", "output_data", "Rabbit_demography_values.txt")
  # demograph_mod            = file.path("rabbit", "output_data", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe")
  # temperature_files        = file.path("data", "environmental_data", "CDS", "kriged_temperature")
  # is_Kelvin                = T
  # precipitation_files      = file.path("data", "environmental_data", "CDS", "precipitation.nc")
  # ndvi_files               = list.files(file.path("data", "environmental_data", "NDVI"), full.names = T)[c(1:3)]
  # flood_files              = file.path("data", "environmental_data", "LAST_monthly_floodmaps")
  # height_files             = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif")
  # out_dir                  = file.path("rabbit", "results")
  # N_sim_samples            = 50
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                      Loading libraries')
  cat('\n-----------------------------------------------------------------')
  
  library(terra)
  library(ncdf4)
  library(dplyr)
  library(EasyABC)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(stats) 
  library(sf)
  library(readr)
  library(coda)
  library(tidyterra)
  library(ggspatial)
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                      Load R helper function')
  cat('\n-----------------------------------------------------------------')
  
  sapply(list.files(file.path("rabbit", "functions"), pattern = ".R", full.names = T), source)
  
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                      Creating output directories')
  cat('\n-----------------------------------------------------------------')
  
  ## Location of files/directory
  out_dir_clim <- file.path(out_dir, "clim_variables")
  
  if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
  if(!dir.exists(out_dir_clim)){dir.create(out_dir_clim, recursive = T)}
  
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                   Set up observational data')
  cat('\n-----------------------------------------------------------------')
  
  obs_df <- readr::read_csv(obs_df) %>%
    pivot_longer(
      cols = c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea",
               "RBD-este", "Puntal", "Marismillas", "Abalario", "Hinojos"),
      names_to  = "transect",
      values_to = "KAI"
    ) %>% rename( Fecha = `Fecha...1`,
                  Long_date = `Fecha...2`) 
  
  if(is.null(years_range)){
    end_year = max(year(as.Date(obs_df$Fecha, format = "%d/%m/%Y")))
    start_year = min(year(as.Date(obs_df$Fecha, format = "%d/%m/%Y")))
  } else {
    end_year = years_range[2]
    start_year = years_range[1]
  }
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                   Load environmental rasters')
  cat('\n-----------------------------------------------------------------')
  
  DT_template = rast(template)
  
  r_temp   = load_and_format_rasters(temperature_files, DT_template)
  r_precip = load_and_format_rasters(precipitation_files, DT_template)
  r_mean_ndvi   = load_and_format_rasters(ndvi_files, DT_template, aggregate = T, fun = "mean")
  r_sd_ndvi     = load_and_format_rasters(ndvi_files, DT_template, aggregate = T, fun = "sd")
  r_flood       = load_and_format_rasters(flood_files, DT_template)
  r_height      = load_and_format_rasters(height_files, DT_template, aggregate = T, fun = "mean")
  
  check_time_span_list(raster_list = list("temperature" = r_temp, 
                                          "precipitation" = r_precip, 
                                          "meanNDVI" = r_mean_ndvi, 
                                          "sdNDVI" = r_sd_ndvi, 
                                          "Fooding" = r_flood), 
                       start_year, end_year)
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                 Format rabbit climate variables')
  cat('\n-----------------------------------------------------------------')
  
  calculate_BM_cDM(temp = r_temp,
                   K = is_Kelvin,
                   precip = r_precip,
                   year_min = start_year,
                   year_max = end_year,
                   result_dir = out_dir_clim,
                   coord_rast = DT_template)
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                 Setting up calibration data')
  cat('\n-----------------------------------------------------------------')
  
  ## Create a named list with the predictors/covariates
  Cov_list = list(
    "ndvi_mean"   = r_mean_ndvi,
    "ndvi_sd"     = r_sd_ndvi,
    "flood"       = r_flood,
    "height"      = r_height
  )
  
  # ----------------------------------------------------------------------------------------
  # 5. Calibration settings
  # ----------------------------------------------------------------------------------------
  
  prior_list <- list(
    c("unif", 1,   10),     # K_max
    c("unif", 0, 1),         # ndvi_optimum
    c("unif", 0, 0.3),       # ndvi_sigma
    c("unif", -3, 3),        # ndvi_beta
    c("unif", 0, 1),         # ndvi_sd_optimum
    c("unif", 0, 0.3),       # ndvi_sd_sigma
    c("unif", -3, 3),        # ndvi_sd_beta
    c("unif", 0, 2),         # heigh_optimum
    c("unif", 0, 1),         # height_sigma
    c("unif", -3, 3),        # height_beta
    c("unif", 1, 5),         # dens_opt
    c("unif", 0.2, 0.5),     # R_lambda
    c("unif", 2,  4)#,       # R_sigma
    # c("normal", 0, 2),     # obs_p
    # c("normal", 0, 2)      # obs_p_ndvi
  )
  
  
  # ----------------------------------------------------------------------------------------
  # 6. Run calibration
  # ----------------------------------------------------------------------------------------
  
  make_model_fn <- function(start_year, end_year) {
    function(par) {
      model_fn(params = par, start_yr = start_year, end_yr = end_year)
    }
  }
  
  # Create the wrapped version with your specific constants
  model_fn_wrapped <- make_model_fn(start_year = 2006, end_year = 2024)
  
  abc_result <- ABC_sequential(
    method              = "Lenormand",
    model               = model_fn_wrapped,
    prior               = prior_list,
    p_acc_min           = 0.20,       # stops when acceptance rate drops below 5%
    summary_stat_target = c(0.20, 0), 
    nb_simul            = 30,         # more manageable
    use_seed            = TRUE,
    verbose             = TRUE,
    n_cluster           = 1
  )
  
  
  # ----------------------------------------------------------------------------------------
  # 7. Save result
  # ----------------------------------------------------------------------------------------
  
  param_df <- as.data.frame(cbind(abc_result$param , abc_result$weights))
  colnames(param_df) <- c("K_max", "ndvi_opt", "ndvi_sigma", "height_opt", "height_sigma",
                          "dens_opt", "R_lambda",
                          "R_sigma", "weights")
  # colnames(param_df) <- c("alpha", "beta_param", "K_max", "dens_opt", "R_lambda", 
  #                         "R_sigma", "weights")
  
  write.csv(param_df, file = file.path(out_dir, "posterior_samples_height_no_beta.csv"))
  saveRDS(abc_result, file = file.path(out_dir, "abc_results_height_no_beta.rds"))
  
  # ----------------------------------------------------------------------------------------
  # 8. Evaluate performance
  # ----------------------------------------------------------------------------------------
  
  # if(!exists("param_df")){param_df <- read.csv(file.path(out_dir, "posterior_samples_height_no_beta.csv"))}
  
  ## To-Do: Update the function to assign the coefficients in the dataframe in function:
  ## "rabbit/functions/compare_sim_to_obs.R L17-24
  historic_sim <- run_simulation(posterior_df  = param_df,
                                 start_year    = start_year,
                                 end_year      = end_year,
                                 N_samples     = N_sim_samples,  # Number of samples
                                 compare_obs   = T)              # Plot comparison of simulated KAI with obs KAI
  
  saveRDS(historic_sim, file = file.path(out_dir, "historic_in_sample_simulation.rds"))
  
  plot(historic_sim$pop_maps$Rabbit_Population_distribution_2024_9$mean)
  plot(historic_sim$pop_maps$Rabbit_Population_distribution_2024_1$mean)
  
  plot_df <- obs_df %>%
    mutate(Fecha = as.Date(Fecha, format = "%d/%m/%Y")) %>%
    select(Fecha, transect, KAI) %>% 
    rename(obs_KAI = KAI) %>%
    right_join(., historic_sim$sim_df)
  
  plot <- ggplot(plot_df) +
    geom_line(aes(x = Fecha, y = KAI, group = id)) +
    geom_line(aes(x = Fecha, y = obs_KAI), colour = "red") +
    facet_wrap(vars(transect))
  
  ggsave(plot, filename = file.path(out_dir, "observation_vs_simulated.png"), width = 5, height = 4)
  
}



cat('\n-----------------------------------------------------------------')
cat('\n                      Loading libraries')
cat('\n-----------------------------------------------------------------')

library(terra)
library(dplyr)
library(EasyABC)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stats)
library(sf)
library(readr)
library(coda)
library(tidyterra)
library(ggspatial)

cat('\n-----------------------------------------------------------------')
cat('\n                      Setting file paths')
cat('\n-----------------------------------------------------------------')

## Location of files/directory
template <- file.path("data", "donana_DT", "environmental_data", "template_raster_500.tif")

obs_transects <- file.path("data", "donana_DT", "rabbit", "Transect_oryctolagus.kml")
obs_df <- file.path("data", "donana_DT", "rabbit", "KAI_Rabbit_Night_2024_v1.csv")

rabbit_demographic_rates <- file.path("data", "donana_DT", "rabbit", "Rabbit_demography_values.txt")
demograph_mod <- file.path("data", "donana_DT", "rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe")

temperature_files_raw <- list.files(file.path("data", "donana_DT", "environmental_data", "CDS"),
                                    pattern = "Krigged_Kriged\\.nc$", full.names = TRUE)
extract_year_month <- function(path) {
  bn <- basename(path)
  yr <- as.integer(sub("^(\\d{4})_(\\d{1,2}).*$", "\\1", bn))
  mo <- as.integer(sub("^(\\d{4})_(\\d{1,2}).*$", "\\2", bn))
  c(yr, mo)
}
if (length(temperature_files_raw) > 0) {
  order_idx <- order(sapply(temperature_files_raw, function(x) extract_year_month(x)[1]),
                     sapply(temperature_files_raw, function(x) extract_year_month(x)[2]))
  temperature_files <- temperature_files_raw[order_idx]
} else {
  temperature_files <- character(0)
}
precipitation_file <- file.path("data", "donana_DT", "environmental_data", "CDS", "precipitation.grib")

ndvi_map   <- file.path("data", "donana_DT", "environmental_data", "NDVI")
flood_map  <- file.path("data", "donana_DT", "environmental_data", "LAST_flooding")

height_high_res <- file.path("data", "donana_DT", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif")
height_file <- file.path("data", "donana_DT", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated_500res.tif")

out_dir_clim <- file.path("data", "donana_DT", "rabbit", "clim_variables")
if (!dir.exists(out_dir_clim)) dir.create(out_dir_clim, recursive = TRUE)

cat('\n-----------------------------------------------------------------')
cat('\n                      Load rabbit helper functions')
cat('\n-----------------------------------------------------------------')

helper_dir <- file.path("rabbit", "functions")
if (!dir.exists(helper_dir)) {
  stop("Helper directory not found: ", helper_dir)
}

sapply(list.files(helper_dir, pattern = "\\.R$", full.names = TRUE), source)

cat('\n-----------------------------------------------------------------')
cat('\n                   Set up observational data')
cat('\n-----------------------------------------------------------------')

obs_df <- readr::read_csv(obs_df) %>%
  pivot_longer(
    cols = c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea",
             "RBD-este", "Puntal", "Marismillas", "Abalario", "Hinojos"),
    names_to  = "transect",
    values_to = "KAI"
  ) %>% rename( Fecha = `Fecha...1`,
                Long_date = `Fecha...2`)

end_year <- max(year(as.Date(obs_df$Fecha, format = "%d/%m/%Y")))
cat('start year hard coded at 2006 for now\n')
start_year <- 2006

DT_template <- rast(template)

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

if (!file.exists(height_file)) {
  height_rast <- rast(height_high_res)
  height_rast <- resample(height_rast, DT_template)
  writeRaster(height_rast, filename = height_file, overwrite = TRUE)
}

cat('\n-----------------------------------------------------------------')
cat('\n                 Format rabbit climate variables')
cat('\n-----------------------------------------------------------------')

calculate_BM_cDM(temp_files = temperature_files,
                 precip_files = precipitation_file,
                 year_min = start_year,
                 year_max = end_year,
                 result_dir = out_dir_clim,
                 coord_rast = DT_template)

cat('\n-----------------------------------------------------------------')
cat('\n         Load posterior calibration results for forecast')
cat('\n-----------------------------------------------------------------')

param_file <- file.path("rabbit", "results", "posterior_samples_height_no_beta.csv")
if (!exists("param_df")) {
  if (file.exists(param_file)) {
    param_df <- read_csv(param_file)
  } else {
    stop("Calibration results file not found: ", param_file)
  }
}

# ----------------------------------------------------------------------------------------
# 9. Creating forecasts for 2027 scenarios
# ----------------------------------------------------------------------------------------

forecast_end_year <- 2027

scenario_definitions <- list(
  warm_dry = list(
    temp_scen = file.path("data", "donana_DT", "ndvi", "results", "temp_warm_values.tif"),
    temp_scen_time = file.path("data", "donana_DT", "ndvi", "results", "temp_warm_time.rds"),
    precip_scen = file.path("data", "donana_DT", "ndvi", "results", "precip_dry_values.tif"),
    precip_scen_time = file.path("data", "donana_DT", "ndvi", "results", "precip_dry_time.rds"),
    ndvi = file.path("data", "donana_DT", "ndvi", "results", "predictions", "NDVI_preset_warm_dry.tif")
  ),
  average = list(
    temp_scen = file.path("data", "donana_DT", "ndvi", "results", "temp_mean_values.tif"),
    temp_scen_time = file.path("data", "donana_DT", "ndvi", "results", "temp_mean_time.rds"),
    precip_scen = file.path("data", "donana_DT", "ndvi", "results", "precip_mean_values.tif"),
    precip_scen_time = file.path("data", "donana_DT", "ndvi", "results", "precip_mean_time.rds"),
    ndvi = file.path("data", "donana_DT", "ndvi", "results", "predictions", "NDVI_preset_all_avg.tif")
  ),
  cold_wet = list(
    temp_scen = file.path("data", "donana_DT", "ndvi", "results", "temp_cold_values.tif"),
    temp_scen_time = file.path("data", "donana_DT", "ndvi", "results", "temp_cold_time.rds"),
    precip_scen = file.path("data", "donana_DT", "ndvi", "results", "precip_wet_values.tif"),
    precip_scen_time = file.path("data", "donana_DT", "ndvi", "results", "precip_wet_time.rds"),
    ndvi = file.path("data", "donana_DT", "ndvi", "results", "predictions", "NDVI_preset_cold_wet.tif")
  )
)

future_scenarios_dir <- file.path("rabbit", "results",)
dir.create(future_scenarios_dir, recursive = TRUE, showWarnings = FALSE)

results_dir <- file.path("data", "donana_DT", "rabbit", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

combined_sim_df <- tibble()

for (variant in variants) {
  for (scenario in names(scenario_definitions)) {
    scenario_info <- scenario_definitions[[scenario]]
    scenario_clim_out <- file.path(future_scenarios_dir, variant, scenario)
    dir.create(scenario_clim_out, recursive = TRUE, showWarnings = FALSE)
    
    climate_files <- build_variant_climate_stacks(
      temp_obs_files = temperature_files,
      precip_obs_file = precipitation_file,
      temp_scen_file = scenario_info$temp_scen,
      temp_scen_time = scenario_info$temp_scen_time,
      precip_scen_file = scenario_info$precip_scen,
      precip_scen_time = scenario_info$precip_scen_time,
      variant = variant,
      scenario = scenario,
      output_dir = scenario_clim_out,
      end_year = forecast_end_year
    )
    
    scenario_clim_dir <- file.path(out_dir_clim, variant, scenario)
    dir.create(scenario_clim_dir, recursive = TRUE, showWarnings = FALSE)
    
    calculate_BM_cDM(
      temp_files = climate_files$temp,
      precip_files = climate_files$precip,
      year_min = start_year,
      year_max = forecast_end_year,
      result_dir = scenario_clim_dir,
      coord_rast = DT_template,
      K = F
    )
    
    pred_ndvi <- rast(scenario_info$ndvi)
    obs_ndvi <- rast(list.files(ndvi_map, full.names = TRUE))
    
    if (variant == "observed_until_2024") {
      ndvi_all <- c(obs_ndvi[[year(time(obs_ndvi)) < 2025]],
                    pred_ndvi[[year(time(pred_ndvi)) >= 2025]])
    } else {
      last_ndvi_date <- max(time(obs_ndvi))
      ndvi_all <- c(obs_ndvi,
                    pred_ndvi[[time(pred_ndvi) > last_ndvi_date]])
    }
    
    ndvi_scene_file <- file.path(results_dir, "ndvi", paste0(variant, "_", scenario, ".tif"))
    if (!dir.exists(dirname(ndvi_scene_file))) dir.create(dirname(ndvi_scene_file), recursive = TRUE)
    writeRaster(ndvi_all, filename = ndvi_scene_file, overwrite = TRUE)
    
    Cov_list <- list(
      ndvi_map = ndvi_scene_file,
      flood_map = flood_map,
      height_file = height_file
    )
    
    sim_output <- run_simulation(
      posterior_df = param_df,
      start_year = start_year,
      end_year = forecast_end_year,
      N_samples = 5,
      compare_obs = TRUE,
      demographic_climate = scenario_clim_dir
    )
    
    pop_out_dir <- file.path(results_dir, "pop_maps", variant, scenario)
    save_population_mean_rasters(sim_output$pop_maps, pop_out_dir, variant, scenario)
    
    combined_sim_df <- bind_rows(
      combined_sim_df,
      sim_output$sim_df %>% mutate(variant = variant, scenario = scenario)
    )
  }
}

write_csv(combined_sim_df, file.path(results_dir, "combined_2027_scenarios_output.csv"))
