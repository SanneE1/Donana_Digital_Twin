#
# Backcast: rabbit densities from the calibrated model on observed environmental data
#
# Separate from rabbit/Workflow_calibration.R and rabbit/Workflow_forecast.R, which are left as
# they are. Uses the posterior of the calibration as is (no recalibration) and the data prepared by
# backcast/Workflow_backcast_data.R (+ backcast/harmonise_landsat_to_s2.R for phase 2).
#
#   phase 1: rabbit_backcast(start_year = 2015)                              Sentinel-2 NDVI
#   phase 2: rabbit_backcast(start_year = 1985,
#                            ndvi_mean_file = "data/backcast/monthly_ndvi_harmonised/mean.tif",
#                            ndvi_sd_file   = "data/backcast/monthly_ndvi_harmonised/sd.tif",
#                            rhd_settings   = modifyList(rhd_default_settings(),
#                                               list(RHD1_arrival_year = 1989, RHD2_arrival_year = 2011)))
#
# rhd_settings: RHD settings for the Pascal program (rabbit/functions/rhd_settings.R), e.g. the
# calibrated ones; NULL = the program's default (schedule mode). Phase 1 starts after both
# arrivals (arrival year 0 in the defaults: RHD present, part of the rabbits immune).
#
# Inputs are the corrected ones (temperature labelled by its own month, gap-filled NDVI from
# Planetary Computer). NB: the posterior was calibrated on the earlier inputs, so until the model
# is recalibrated on these, parameters and inputs do not quite match.
#
# The first simulated year is spin-up (the model starts at carrying capacity) and is left out of
# the summaries. Census dates before the calibration window were not used to calibrate, so they
# test the backcast out of sample.

rabbit_backcast <- function(
    start_year,
    end_year                 = lubridate::year(Sys.Date()) - 1,
    template                 = file.path("data", "environmental_data", "template_raster_500.tif"),

    obs_transects            = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml"),
    obs_df                   = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"),

    rabbit_demographic_rates = file.path("data", "rabbit", "Rabbit_demography_values.txt"),
    demograph_mod            = file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"),
    posterior_file           = file.path("rabbit", "results", "posterior_samples.csv"),
    calibration_years        = c(2017, 2024),

    temperature_files        = file.path("data", "environmental_data", "temperature_CDS"),
    is_Kelvin                = T,
    precipitation_files      = file.path("data", "environmental_data", "precipitation_CDS"),
    ndvi_mean_file           = file.path("data", "environmental_data", "NDVI_large", "harmonised", "mean.tif"),
    ndvi_sd_file             = file.path("data", "environmental_data", "NDVI_large", "harmonised", "sd.tif"),
    flood_files              = file.path("data", "backcast", "monthly_floodmaps"),
    height_files             = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif"),

    out_dir                  = file.path("rabbit", "results", paste0("backcast_", start_year)),
    N_sim_samples            = 50,
    rhd_settings             = NULL
){

  cat('\n-----------------------------------------------------------------')
  cat('\n                      Loading libraries')
  cat('\n-----------------------------------------------------------------')

  library(terra)
  library(ncdf4)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(sf)
  library(readr)

  cat('\n-----------------------------------------------------------------')
  cat('\n                      Load R helper function')
  cat('\n-----------------------------------------------------------------')

  sapply(list.files(file.path("rabbit", "functions"), pattern = "\\.R$", full.names = T), source)
  source(file.path("backcast", "functions", "backcast_helpers.R"))

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

  cat('\n-----------------------------------------------------------------')
  cat('\n                   Load environmental rasters')
  cat('\n-----------------------------------------------------------------')

  DT_template <- rast(template)

  # Climate: ERA5-Land monthly means (data/environmental_data/<variable>_CDS)
  r_temp   <- load_monthly_stack(temperature_files, DT_template, pattern = "\\.nc$")
  r_precip <- load_monthly_stack(precipitation_files, DT_template, pattern = "\\.nc$")

  # NDVI is already aggregated to the template and gap-filled; months before the first image (only
  # the spin-up year of phase 1) get the average of the same month in the first years
  start_date  <- as.Date(sprintf("%d-01-01", start_year))
  r_mean_ndvi <- fill_leading_climatology(read_monthly_ndvi(ndvi_mean_file, DT_template), start_date)
  r_sd_ndvi   <- fill_leading_climatology(read_monthly_ndvi(ndvi_sd_file, DT_template), start_date)

  r_flood  <- load_and_format_rasters(flood_files, DT_template,
                                      aggregate = T, fun = "max",
                                      temporal_agg = "month", temporal_fun = "max")
  r_height <- load_and_format_rasters(height_files, DT_template, aggregate = T, fun = "mean")

  # Nov/Dec before the first year start the dry-month count in calculate_BM_cDM()
  window <- function(r){
    t <- as.Date(time(r))
    r[[t >= as.Date(sprintf("%d-11-01", start_year - 1)) & t <= as.Date(sprintf("%d-12-01", end_year))]]
  }
  r_temp <- window(r_temp); r_precip <- window(r_precip)
  r_mean_ndvi <- window(r_mean_ndvi); r_sd_ndvi <- window(r_sd_ndvi); r_flood <- window(r_flood)

  for(clim in list(r_temp, r_precip)){
    if(min(as.Date(time(clim))) > as.Date(sprintf("%d-11-01", start_year - 1))){
      stop("Climate data must start in November ", start_year - 1, " (found ", min(as.Date(time(clim))), ")")
    }
  }

  # Missing climate months are filled by temporal interpolation (NDVI arrives gap-filled; any cell
  # still without NDVI stops the run in create_K_map_gaussian()). Months without a flood map are
  # treated as not flooded by the K maps, as in the calibration.
  checked_rasters <- check_time_span_list(raster_list = list("temperature"   = r_temp,
                                                             "precipitation" = r_precip,
                                                             "meanNDVI"      = r_mean_ndvi,
                                                             "sdNDVI"        = r_sd_ndvi),
                                          start_year, end_year,
                                          fill_gaps = TRUE)
  r_temp      <- checked_rasters$temperature
  r_precip    <- checked_rasters$precipitation
  r_mean_ndvi <- checked_rasters$meanNDVI
  r_sd_ndvi   <- checked_rasters$sdNDVI

  flood_years <- unique(year(as.Date(time(r_flood))))
  no_flood_years <- setdiff(start_year:end_year, flood_years)
  if(length(no_flood_years) > 0){
    warning("No flood maps at all for: ", paste(no_flood_years, collapse = ", "), call. = FALSE)
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

  Cov_list = list(
    "ndvi_mean"   = r_mean_ndvi,
    "ndvi_sd"     = r_sd_ndvi,
    "flood"       = r_flood,
    "height"      = r_height
  )

  cat('\n-----------------------------------------------------------------')
  cat('\n                 Running the backcast')
  cat('\n-----------------------------------------------------------------')

  param_df <- read.csv(posterior_file)

  backcast_sim <- run_simulation(posterior_df           = param_df,
                                 start_year             = start_year,
                                 end_year               = end_year,
                                 N_samples              = N_sim_samples,
                                 compare_obs            = T,
                                 Cov_list               = Cov_list,
                                 DT_template            = DT_template,
                                 obs_transects          = obs_transects,
                                 obs_df                 = obs_df,
                                 demograph_mod          = demograph_mod,
                                 rabbit_demographic_rates = rabbit_demographic_rates,
                                 demographic_climate    = normalizePath(out_dir_clim),
                                 rhd_settings           = rhd_settings)

  save_simulation(backcast_sim, file.path(out_dir, "simulation"))
  if (!is.null(backcast_sim$rhd)) {
    write_csv(backcast_sim$rhd$outbreaks, file.path(out_dir, "rhd_outbreaks.csv"))
    write_csv(backcast_sim$rhd$monthly, file.path(out_dir, "rhd_monthly.csv"))
  }

  cat('\n-----------------------------------------------------------------')
  cat('\n                 Summaries')
  cat('\n-----------------------------------------------------------------')

  spin_up_end <- as.Date(sprintf("%d-12-31", start_year))
  calib_start <- as.Date(sprintf("%d-01-01", calibration_years[1]))
  calib_end   <- as.Date(sprintf("%d-12-31", calibration_years[2]))

  # Total abundance over the study area per mapped month. The interval is the sum of the
  # cell-wise 95% intervals, so wider than the interval of the total itself.
  total_df <- bind_rows(lapply(names(backcast_sim$pop_maps), function(nm){
    ym <- as.integer(regmatches(nm, regexpr("\\d{4}_\\d{1,2}$", nm)) %>% strsplit("_") %>% unlist())
    m  <- backcast_sim$pop_maps[[nm]]
    data.frame(date     = as.Date(sprintf("%d-%02d-15", ym[1], ym[2])),
               total    = global(m$mean, "sum", na.rm = TRUE)[, 1],
               lower_ci = global(m$lower_ci, "sum", na.rm = TRUE)[, 1],
               upper_ci = global(m$upper_ci, "sum", na.rm = TRUE)[, 1])
  })) %>%
    filter(date > spin_up_end) %>%
    arrange(date)

  write_csv(total_df, file.path(out_dir, "total_abundance.csv"))

  total_plot <- ggplot(total_df, aes(date, total)) +
    annotate("rect", xmin = calib_start, xmax = calib_end, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#c87941", alpha = 0.25) +
    geom_line(colour = "#8B4513") +
    geom_point(colour = "#8B4513", size = 0.8) +
    labs(x = NULL, y = "Simulated rabbits, study area total",
         title = paste("Rabbit backcast", start_year + 1, "to", end_year),
         subtitle = "Grey: calibration period. Band: summed cell-wise 95% intervals.") +
    theme_minimal(base_size = 12)
  ggsave(total_plot, filename = file.path(out_dir, "total_abundance.png"), width = 8, height = 4)

  # Simulated vs observed KAI per transect. As in the calibration, the observed KAI is rescaled
  # per transect, with scaling rates from the calibration period only.
  obs_plot_df <- obs_df %>%
    mutate(Fecha = as.Date(Fecha, format = "%d/%m/%Y")) %>%
    select(Fecha, transect, KAI) %>%
    filter(!is.na(KAI), Fecha > spin_up_end)

  sim_df <- backcast_sim$sim_df %>% filter(Fecha > spin_up_end)

  sim_sum <- sim_df %>%
    group_by(transect, Fecha) %>%
    summarise(mean = mean(KAI), min = min(KAI), max = max(KAI), .groups = "drop")

  rates <- sim_sum %>%
    inner_join(obs_plot_df, by = c("Fecha", "transect")) %>%
    filter(KAI > 0, Fecha >= calib_start, Fecha <= calib_end) %>%
    group_by(transect) %>%
    summarise(rate = mean(mean / KAI, na.rm = TRUE), .groups = "drop")

  obs_scaled <- obs_plot_df %>% left_join(rates, by = "transect") %>% mutate(KAI = KAI * rate)

  kai_plot <- ggplot() +
    annotate("rect", xmin = calib_start, xmax = calib_end, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_ribbon(data = sim_sum, aes(x = Fecha, ymin = min, ymax = max),
                fill = "#c87941", alpha = 0.25) +
    geom_line(data = sim_sum, aes(x = Fecha, y = mean), colour = "#8B4513", linewidth = 0.8) +
    geom_line(data = obs_scaled, aes(x = Fecha, y = KAI), colour = "black", linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(vars(transect), scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Relative rabbit abundance (KAI)",
         title = "Backcast vs census per transect",
         subtitle = paste("Brown: simulated mean and range of", N_sim_samples, "posterior draws;",
                          "black: observed KAI (rescaled); grey: calibration period")) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"),
          strip.text = element_text(face = "bold"))
  ggsave(kai_plot, filename = file.path(out_dir, "backcast_vs_census.png"), width = 9, height = 9)

  # Agreement in and out of the calibration period
  agreement <- sim_sum %>%
    inner_join(obs_plot_df, by = c("Fecha", "transect")) %>%
    mutate(period = ifelse(Fecha >= calib_start & Fecha <= calib_end, "calibration", "out of sample")) %>%
    group_by(transect, period) %>%
    summarise(n_censuses = n(),
              spearman = suppressWarnings(cor(mean, KAI, method = "spearman")),
              .groups = "drop")
  write_csv(agreement, file.path(out_dir, "agreement_with_census.csv"))
  print(agreement)

  invisible(list(sim = backcast_sim, total = total_df, agreement = agreement))
}
