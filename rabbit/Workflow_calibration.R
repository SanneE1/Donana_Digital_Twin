# Rabbit forecasting calibration

rabbit_calibration <- function(
  template       = file.path("data", "environmental_data", "template_raster_500.tif"),
  years_range    = NULL,

  obs_transects  = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml"),
  obs_df         = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"),

  rabbit_demographic_rates = file.path("rabbit", "output_data", "Rabbit_demography_values.txt"),
  demograph_mod            = file.path("rabbit", "output_data", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"),

  temperature_files   = file.path("data", "environmental_data", "CDS", "kriged_temperature"),
  is_Kelvin           = T,
  precipitation_files = file.path("data", "environmental_data", "CDS", "precipitation.nc"),
  ndvi_files          = list.files(file.path("data", "environmental_data", "NDVI"), full.names = T),
  flood_files         = file.path("data", "environmental_data", "LAST_monthly_floodmaps"),
  height_files        = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif"),

  out_dir             = file.path("rabbit", "results"),
  N_sim_samples       = 50

){

  # template                 = file.path("data", "environmental_data", "template_raster_500.tif")
  # years_range              = c(2017,2024)
  # obs_transects            = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml")
  # obs_df                   = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv")
  # rabbit_demographic_rates = file.path("rabbit", "output_data", "Rabbit_demography_values.txt")
  # demograph_mod            = file.path("rabbit", "output_data", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe")
  # temperature_files        = file.path("data", "environmental_data", "CDS", "kriged_temperature")
  # is_Kelvin           = T
  # precipitation_files = file.path("data", "environmental_data", "CDS", "precipitation.nc")
  # # ndvi_files               = file.path("data", "environmental_data", "NDVI")
  # ndvi_files               = list.files(file.path("data", "environmental_data", "NDVI"), full.names = T)
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
  
  r_temp   = load_and_format_rasters(temperature_files, DT_template, 
                                      temporal_agg = "month", temporal_fun = "mean")
  r_precip = load_and_format_rasters(precipitation_files, DT_template, temporal_agg = "month", 
                                      temporal_fun = "sum")
  r_mean_ndvi   = load_and_format_rasters(ndvi_files, DT_template, 
                                          aggregate = T, fun = "mean",
                                          temporal_agg = "month", temporal_fun = "max")
  r_sd_ndvi     = load_and_format_rasters(ndvi_files, DT_template, 
                                          aggregate = T, fun = "sd",
                                          temporal_agg = "month", temporal_fun = "max")
  r_flood       = load_and_format_rasters(flood_files, DT_template,
                                          aggregate = T, fun = "max",
                                          temporal_agg = "month", temporal_fun = "max")
  r_height      = load_and_format_rasters(height_files, DT_template, aggregate = T, fun = "mean")
  
  check_time_span_list(raster_list = list("temperature" = r_temp, 
                                          "precipitation" = r_precip, 
                                          "meanNDVI" = r_mean_ndvi, 
                                          "sdNDVI" = r_sd_ndvi, 
                                          "Fooding" = r_flood #, 
                                          #"VegHeight" = r_height
                                          ), 
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
    c("unif", 1, 30),        # K_max
    c("unif", 0, 1),         # ndvi_optimum
    c("unif", 0, 0.3),       # ndvi_sigma
    c("unif", 0, 1),         # ndvi_sd_optimum
    c("unif", 0, 0.3),       # ndvi_sd_sigma
    c("unif", 0, 2),         # heigh_optimum
    c("unif", 0, 1),         # height_sigma
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
  model_fn_wrapped <- make_model_fn(start_year = start_year, end_year = end_year)
  
  abc_result <- ABC_sequential(
    method              = "Lenormand",
    model               = model_fn_wrapped,
    prior               = prior_list,
    p_acc_min           = 0.20,       # stops when acceptance rate drops below 5%
    summary_stat_target = c(0.20, 0), 
    nb_simul            = 100,         
    use_seed            = TRUE,
    verbose             = TRUE,
    n_cluster           = 1
  )
  
  
  # ----------------------------------------------------------------------------------------
  # 7. Save result
  # ----------------------------------------------------------------------------------------
  
  param_df <- as.data.frame(cbind(abc_result$param , abc_result$weights))
  colnames(param_df) <- c("K_max", "ndvi_opt", "ndvi_sigma",
                          "ndvi_sd_opt", "ndvi_sd_sigma",
                          "height_opt", "height_sigma",
                          "dens_opt", "R_lambda", "R_sigma", "weights")
  
  write.csv(param_df, file = file.path(out_dir, "posterior_samples_height_no_beta.csv"))
  saveRDS(abc_result, file = file.path(out_dir, "abc_results_height_no_beta.rds"))
  
  # ----------------------------------------------------------------------------------------
  # 8. Evaluate performance
  # ----------------------------------------------------------------------------------------
  
  if(!exists("param_df")){param_df <- read.csv(file.path(out_dir, "posterior_samples_height_no_beta.csv"))}
  
  historic_sim <- run_simulation(posterior_df           = param_df,
                                 start_year             = start_year,
                                 end_year               = end_year,
                                 N_samples              = N_sim_samples,  # Number of samples
                                 compare_obs            = T,             # Plot comparison of simulated KAI with obs KAI
                                 Cov_list               = Cov_list,
                                 DT_template            = DT_template,
                                 obs_transects          = obs_transects,
                                 obs_df                 = obs_df,
                                 demograph_mod          = demograph_mod,
                                 rabbit_demographic_rates = rabbit_demographic_rates,
                                 demographic_climate    = out_dir_clim)
  
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




rabbit_calibration(use_years = c(2016, 2024))


