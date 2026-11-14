#
# Rabbit forecasting calibration
#

rabbit_calibration <- function(
    template       = file.path("data", "environmental_data", "template_raster_500.tif"),
    years_range    = c(2019, 2022),
    
    obs_transects  = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml"),
    obs_df         = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"),
    
    rabbit_demographic_rates = file.path("rabbit", "output_data", "Rabbit_demography_values.txt"),
    demograph_mod            = file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"),
    
    temperature_files   = file.path("data", "environmental_data", "temperature_CDS"),
    is_Kelvin           = T,
    precipitation_files = file.path("data", "environmental_data", "precipitation_CDS"),
    ndvi_mean_file          = file.path("data", "environmental_data", "NDVI_large", "harmonised", "mean.tif"),
    ndvi_sd_file            = file.path("data", "environmental_data", "NDVI_large", "harmonised", "sd.tif"),
    flood_files         = file.path("data", "environmental_data", "LAST_monthly_floodmaps"),
    height_files        = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif"),
    
    out_dir             = file.path("rabbit", "results"),
    N_sim_samples       = 50,
    n_cluster           = 4,
    fill_gaps           = TRUE   # if TRUE, fill missing monthly raster layers via temporal interpolation
    
){
  # 
  #   template                 = file.path("data", "environmental_data", "template_raster_500.tif")
  #   years_range              = c(2006,2024)
  #   obs_transects            = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml")
  #   obs_df                   = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv")
  # rabbit_demographic_rates = file.path("data", "rabbit", "Rabbit_demography_values.txt")
  # demograph_mod            = file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe")
  #   temperature_files        = file.path("data", "environmental_data", "temperature_CDS")
  #   is_Kelvin           = T
  #   precipitation_files = file.path("data", "environmental_data", "precipitation_CDS")
  #   # ndvi_files               = file.path("data", "environmental_data", "NDVI")
  #   ndvi_files               = file.path("data", "environmental_data", "monthly_ndvi.tif")
  #   flood_files              = file.path("data", "environmental_data", "LAST_monthly_floodmaps")
  #   height_files             = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif")
  #   out_dir                  = file.path("rabbit", "results")
  #   N_sim_samples            = 50
  #   n_cluster                = 4
  #   fill_gaps                = TRUE
  
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
  library(ggh4x)
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                      Load R helper function')
  cat('\n-----------------------------------------------------------------')
  
  # Captured as absolute paths so parallel workers (which start in a fresh
  # R session, possibly with a different working directory) can still find
  # and source these functions.
  proj_root    <- normalizePath(getwd())
  function_dir <- normalizePath(file.path("rabbit", "functions"))
  
  sapply(list.files(function_dir, pattern = ".R", full.names = T), source)
  
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                      Creating output directories')
  cat('\n-----------------------------------------------------------------')
  
  ## Location of files/directory
  out_dir_clim <- file.path(out_dir, "clim_variables")
  
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
  # NDVI mean and within-cell SD: monthly 500 m products (NDVI_large/harmonised), cut to the template
  r_mean_ndvi   = load_and_format_rasters(ndvi_mean_file, DT_template,
                                          temporal_agg = "month", temporal_fun = "max")
  r_sd_ndvi     = load_and_format_rasters(ndvi_sd_file, DT_template,
                                          temporal_agg = "month", temporal_fun = "max")
  r_flood       = load_and_format_rasters(flood_files, DT_template,
                                          aggregate = T, fun = "max",
                                          temporal_agg = "month", temporal_fun = "max")
  r_height      = load_and_format_rasters(height_files, DT_template, aggregate = T, fun = "mean")
  
  checked_rasters <- check_time_span_list(raster_list = list("temperature" = r_temp,
                                                             "precipitation" = r_precip,
                                                             "meanNDVI" = r_mean_ndvi,
                                                             "sdNDVI" = r_sd_ndvi,
                                                             "Fooding" = r_flood #,
                                                             #"VegHeight" = r_height
  ),
  start_year, end_year,
  fill_gaps = fill_gaps)
  
  if(fill_gaps){
    r_temp      <- checked_rasters$temperature
    r_precip    <- checked_rasters$precipitation
    r_mean_ndvi <- checked_rasters$meanNDVI
    r_sd_ndvi   <- checked_rasters$sdNDVI
    #r_flood     <- checked_rasters$Fooding
  }
  
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
    c("unif", 15, 45),       # K_max
    c("unif", 0.1, 0.9),     # ndvi_optimum
    c("unif", 0.01, 0.5),    # ndvi_sigma
    c("unif", 0.1, 0.9),     # ndvi_sd_optimum
    c("unif", 0.01, 0.5),    # ndvi_sd_sigma
    c("unif", 0, 2),         # heigh_optimum
    c("unif", 0.01, 1),      # height_sigma
    c("unif", 3, 20),        # dens_opt
    c("unif", 0.2, 0.5),     # R_lambda
    c("unif", 2,  4),        # R_sigma
    c("unif", 0, 1),         # ndvi_w (weights are normalised to sum to 1)
    c("unif", 0, 1),         # ndvi_sd_w
    c("unif", 0, 1)#,        # height_w
    # c("normal", 0, 2),     # obs_p
    # c("normal", 0, 2)      # obs_p_ndvi
  )
  
  # RHD, immunity mode (see docs/model_settings.md and rabbit/functions/rhd_settings.R). Both strains
  # arrived before the calibration period, so the disease is present from the start (arrival year 0
  # in the defaults), with part of the rabbits immune. Settings not listed in rhd_priors keep the
  # values of rhd_default_settings().
  rhd_settings <- rhd_default_settings()
  
  # Calibrated RHD settings; their values follow the 13 parameters above (params[15], ...).
  # A "log10_" name is calibrated on the log10 scale.
  rhd_priors <- list(
    RHD_f_crit      = c("unif", 0.2, 0.9),    # susceptible share at which an exposed cell has a 50% outbreak chance
    log10_RHD_spark = c("unif", -4, -1.5),    # background introduction, per cell per month (1e-4 to 0.03)
    RHD_spread      = c("unif", 0, 1),        # spread per neighbouring cell with an outbreak last month
    RHD_mortality   = c("unif", 0.3, 0.9)     # share of infected rabbits that die
  )
  prior_list <- c(prior_list, unname(rhd_priors))
  
  # Third summary statistic: how far the median interval between park-wide outbreak waves lies
  # outside 24-36 months (0 = inside; see rhd_interval_distance). The first (spin-up) year is left out.
  rhd_target <- list(scale = "park", target = c(24, 36), from = rhd_month_index(start_year + 1, 1))
  
  
  # ----------------------------------------------------------------------------------------
  # 6. Run calibration
  # ----------------------------------------------------------------------------------------
  
  # Absolute paths so the model still finds its inputs when a parallel
  # worker starts in a working directory other than proj_root.
  demograph_mod_abs            <- normalizePath(demograph_mod)
  check_rabbit_model(demograph_mod_abs)   # needs the build with named arguments
  rabbit_demographic_rates_abs <- normalizePath(rabbit_demographic_rates)
  obs_transects_abs            <- normalizePath(obs_transects)
  out_dir_clim_abs             <- normalizePath(out_dir_clim)
  
  # terra SpatRasters hold external C pointers that cannot be serialized to
  # other R processes as-is; wrap()/unwrap() convert them to/from a portable
  # form so they can be shipped to n_cluster > 1 parallel workers.
  Cov_list_wrapped    <- lapply(Cov_list, terra::wrap)
  DT_template_wrapped <- terra::wrap(DT_template)
  
  make_model_fn <- function(start_year, end_year,
                            Cov_list_wrapped, DT_template_wrapped,
                            demograph_mod, rabbit_demographic_rates,
                            obs_transects, obs_df, out_dir_clim,
                            function_dir, proj_root,
                            rhd_settings, rhd_param_names, rhd_target) {
    # Function arguments are lazy in R: until referenced, each one above is
    # just a promise pointing back at THIS call's expressions, evaluated in
    # rabbit_calibration()'s frame. That frame also holds large objects
    # (unwrapped SpatRasters, etc.) that don't survive being sent to a
    # parallel worker, so if those promises are forced there instead of
    # here, resolving them on the worker fails ("object ... not found").
    # Force everything now, then give the closure its own clean environment
    # containing only these already-resolved values, so nothing from
    # rabbit_calibration()'s frame is dragged along when the closure is
    # serialized to a worker.
    force(start_year); force(end_year)
    force(Cov_list_wrapped); force(DT_template_wrapped)
    force(demograph_mod); force(rabbit_demographic_rates)
    force(obs_transects); force(obs_df); force(out_dir_clim)
    force(function_dir); force(proj_root)
    force(rhd_settings); force(rhd_param_names); force(rhd_target)
    
    model_fn_inner <- function(par) {
      # A parallel (PSOCK) worker starts as a blank R session with no
      # packages loaded, no sourced functions, and no access to this
      # script's environment, so every call has to be self-sufficient.
      setwd(proj_root)
      library(terra); library(dplyr); library(tidyr)
      library(lubridate); library(sf); library(readr)
      if(!exists("model_fn", mode = "function")){
        sapply(list.files(function_dir, pattern = ".R", full.names = T), source)
      }
      
      model_fn(params        = par,
               start_yr      = start_year,
               end_yr        = end_year,
               Cov_list      = lapply(Cov_list_wrapped, terra::unwrap),
               DT_template   = terra::unwrap(DT_template_wrapped),
               demograph_mod = demograph_mod,
               rabbit_demographic_rates = rabbit_demographic_rates,
               obs_transects = obs_transects,
               obs_df        = obs_df,
               out_dir_clim  = out_dir_clim,
               rhd_settings  = rhd_settings,
               rhd_param_names = rhd_param_names,
               rhd_target    = rhd_target)
    }
    
    environment(model_fn_inner) <- list2env(
      list(start_year = start_year, end_year = end_year,
           Cov_list_wrapped = Cov_list_wrapped, DT_template_wrapped = DT_template_wrapped,
           demograph_mod = demograph_mod, rabbit_demographic_rates = rabbit_demographic_rates,
           obs_transects = obs_transects, obs_df = obs_df, out_dir_clim = out_dir_clim,
           function_dir = function_dir, proj_root = proj_root,
           rhd_settings = rhd_settings, rhd_param_names = rhd_param_names, rhd_target = rhd_target),
      parent = globalenv()
    )
    model_fn_inner
  }
  
  # Create the wrapped version with your specific constants
  model_fn_wrapped <- make_model_fn(start_year           = start_year,
                                    end_year              = end_year,
                                    Cov_list_wrapped      = Cov_list_wrapped,
                                    DT_template_wrapped   = DT_template_wrapped,
                                    demograph_mod         = demograph_mod_abs,
                                    rabbit_demographic_rates = rabbit_demographic_rates_abs,
                                    obs_transects         = obs_transects_abs,
                                    obs_df                = obs_df,
                                    out_dir_clim          = out_dir_clim_abs,
                                    function_dir          = function_dir,
                                    proj_root             = proj_root,
                                    rhd_settings          = rhd_settings,
                                    rhd_param_names       = names(rhd_priors),
                                    rhd_target            = rhd_target)
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                 Running ABC')
  cat('\n-----------------------------------------------------------------')
  
  abc_result <- ABC_sequential(
    method              = "Lenormand",
    model               = model_fn_wrapped,
    prior               = prior_list,
    inside_prior        = FALSE,
    p_acc_min           = 0.05,       # stops when acceptance rate drops below 5%
    # model_fn returns: median correlation of simulated and observed KAI over the transects (best: 1),
    # median correlation of their changes between censuses (best: 1), and the RHD interval distance
    # (best: 0). NB: the target used to be c(0, 0), i.e. no correlation with the observations.
    summary_stat_target = c(1, 1, 0),
    nb_simul            = 300,
    use_seed            = TRUE,
    verbose             = TRUE,
    n_cluster           = n_cluster
  )
  
  
  
  
  # ----------------------------------------------------------------------------------------
  # 7. Save result
  # ----------------------------------------------------------------------------------------
  
  cat('\n-----------------------------------------------------------------')
  cat('\n                 Saving results')
  cat('\n-----------------------------------------------------------------')
  
  param_df <- as.data.frame(cbind(abc_result$param , abc_result$weights))
  colnames(param_df) <- c("K_max", "ndvi_opt", "ndvi_sigma",
                          "ndvi_sd_opt", "ndvi_sd_sigma",
                          "height_opt", "height_sigma",
                          "dens_opt", "R_lambda", "R_sigma",
                          "ndvi_w", "ndvi_sd_w", "height_w", names(rhd_priors), "weights")
  
  write.csv(param_df, file = file.path(out_dir, "posterior_samples.csv"))
  saveRDS(abc_result, file = file.path(out_dir, "posterior_samples.rds"))
  
  # ----------------------------------------------------------------------------------------
  # 8. Evaluate performance
  # ----------------------------------------------------------------------------------------
  
  if(!exists("param_df")){param_df <- read.csv(file.path(out_dir, "posterior_samplesc.csv"))}
  
  historic_sim <- run_simulation(posterior_df           = param_df,
                                 start_year             = 2019,
                                 end_year               = 2022,
                                 N_samples              = N_sim_samples,  # Number of samples
                                 compare_obs            = T,              # Plot comparison of simulated KAI with obs KAI
                                 Cov_list               = Cov_list,
                                 DT_template            = DT_template,
                                 obs_transects          = obs_transects,
                                 obs_df                 = obs_df,
                                 demograph_mod          = demograph_mod,
                                 rabbit_demographic_rates = rabbit_demographic_rates,
                                 demographic_climate    = out_dir_clim,
                                 rhd_settings           = rhd_settings)
  
  save_simulation(historic_sim, file.path(out_dir, "historic_in_sample_simulation"))
  historic_sim <- load_simulation(file.path(out_dir, "historic_in_sample_simulation"))
  
  plot(historic_sim$pop_maps$Rabbit_Population_distribution_2022_9$mean)
  plot(historic_sim$pop_maps$Rabbit_Population_distribution_2020_9$mean)
  
  sim_sum <- historic_sim$sim_df %>%
    group_by(transect, Fecha) %>%
    summarise(mean = mean(KAI), min = min(KAI), max = max(KAI), .groups = "drop")
  
  rates <- historic_sim$sim_df %>%
    rename(KAI_sim = KAI) %>%
    inner_join(obs_df %>% mutate(Fecha = as.Date(obs_df$Fecha, format = "%d/%m/%Y")), 
               by = c("Fecha", "transect")) %>%
    filter(KAI > 0) %>%
    mutate(rate = KAI_sim / KAI) %>%
    group_by(transect) %>%
    summarise(rate = mean(rate, na.rm = TRUE), .groups = "drop")
  
  scaled <- function(d) d %>% left_join(rates, by = "transect") %>% mutate(KAI = KAI * rate)
  
  # y-limit per transect: simulation max (with 5% headroom), floored at 0
  lims <- historic_sim$sim_df %>%
    group_by(transect) %>%
    summarise(hi = max(KAI, na.rm = TRUE) * 1.05, .groups = "drop") %>%
    arrange(transect)
  
  y_scales <- lapply(lims$hi, function(h) scale_y_continuous(limits = c(0, h)))
  
  plot <- ggplot() +
    geom_vline(xintercept = as.Date(max(obs_df$Fecha, na.rm = T), format = "%d/%m/%Y"),
               linetype = "dotted", colour = "grey50") +
    geom_ribbon(data = sim_sum, aes(x = Fecha, ymin = min, ymax = max),
                fill = "#c87941", alpha = 0.25) +
    geom_line(data = historic_sim$sim_df, aes(x = Fecha, y = KAI, group = id),
              colour = "#8B4513", alpha = 0.25, linewidth = 0.3) +
    geom_line(data = sim_sum, aes(x = Fecha, y = mean), colour = "#8B4513", linewidth = 1) +
    geom_line(data = scaled(obs_df) %>% mutate(Fecha = as.Date(Fecha, format = "%d/%m/%Y")),
              aes(x = Fecha, y = KAI), colour = "black", linewidth = 0.6, na.rm = TRUE) +
    facet_wrap(vars(transect), scales = "free_y", ncol = 2) +
    facetted_pos_scales(y = y_scales) +
    coord_cartesian(xlim = c(as.Date("2022-01-01"), max(historic_sim$sim_df$Fecha))) +
    labs(x = NULL, y = "Relative rabbit abundance (KAI)",
         title = "Rabbit abundance per transect",
         subtitle = paste("Brown: simulated mean (line), range of 5 samples (band);",
                          "black: observed KAI (rescaled); blue: 2025 census.")) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"),
          strip.text = element_text(face = "bold"))
  
  ggsave(plot, filename = file.path(out_dir, "observation_vs_simulated.png"), width = 5, height = 4)
  
}




rabbit_calibration(years_range = c(2017, 2024))


