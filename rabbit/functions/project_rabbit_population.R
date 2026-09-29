# project_rabbit_population(
#   template        = file.path("data", "environmental_data", "template_raster_500.tif"),
#   obs_df          = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "KAI_Rabbit_Night_2024_v1.csv"),
#   covariate_files = list(
#     temperature_files   = file.path("data", "environmental_data", "CDS", "kriged_temperature"),
#     precipitation_files = file.path("data", "environmental_data", "CDS", "precipitation.nc"),
#     ndvi_files           = list.files(file.path("data", "environmental_data", "NDVI"), full.names = T),
#     flood_files          = file.path("data", "environmental_data", "LAST_monthly_floodmaps"),
#     height_files         = file.path("data", "environmental_data", "vegetation_height", "merged_veg_height_donana_NA_interpolated.tif")
#   ),
#   out_dir = file.path("rabbit", "results")
# )

# Statistical (non-mechanistic) spatial projection of rabbit abundance:
# fits KAI observed at the census transects against habitat covariates and
# the climate of the `climate_window_months` months leading up to each
# census, then predicts that relationship across every cell of the template
# raster for the most recent census date. Writes two GeoTIFFs: the predicted
# mean and the associated prediction uncertainty (standard error) per cell.

project_rabbit_population <- function(
    template,
    obs_df,
    covariate_files,
    out_dir,
    obs_transects          = file.path("data", "rabbit", "Doñana_census_data_PacoCarro", "Transect_oryctolagus.kml"),
    is_Kelvin              = TRUE,
    climate_window_months  = 6,
    mean_filename           = "Rabbit_population_mean.tif",
    uncertainty_filename    = "Rabbit_population_uncertainty.tif"
) {

  required_cov <- c("temperature_files", "precipitation_files", "ndvi_files", "flood_files", "height_files")
  missing_cov  <- setdiff(required_cov, names(covariate_files))
  if (length(missing_cov) > 0) {
    stop("covariate_files is missing: ", paste(missing_cov, collapse = ", "))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat('\n-----------------------------------------------------------------')
  cat('\n         Loading covariate and climate rasters')
  cat('\n-----------------------------------------------------------------\n')

  DT_template <- rast(template)

  r_temp <- load_and_format_rasters(covariate_files$temperature_files, DT_template)
  if (is_Kelvin) r_temp <- r_temp - 273.15
  time(r_temp) <- as.Date(get_stack_time(r_temp))

  r_precip <- load_and_format_rasters(covariate_files$precipitation_files, DT_template)
  time(r_precip) <- as.Date(get_stack_time(r_precip))

  r_ndvi_mean <- load_and_format_rasters(covariate_files$ndvi_files, DT_template, aggregate = TRUE, fun = "mean")
  time(r_ndvi_mean) <- as.Date(get_stack_time(r_ndvi_mean))

  r_ndvi_sd <- load_and_format_rasters(covariate_files$ndvi_files, DT_template, aggregate = TRUE, fun = "sd")
  time(r_ndvi_sd) <- as.Date(get_stack_time(r_ndvi_sd))

  r_flood <- load_and_format_rasters(covariate_files$flood_files, DT_template)
  time(r_flood) <- as.Date(get_stack_time(r_flood))

  r_height <- load_and_format_rasters(covariate_files$height_files, DT_template, aggregate = TRUE, fun = "mean")
  if (nlyr(r_height) > 1) r_height <- r_height[[1]]

  cat('\n-----------------------------------------------------------------')
  cat('\n         Set up observational data')
  cat('\n-----------------------------------------------------------------\n')

  obs <- readr::read_csv(obs_df, show_col_types = FALSE) %>%
    pivot_longer(
      cols = c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea",
               "RBD-este", "Puntal", "Marismillas", "Abalario", "Hinojos"),
      names_to  = "transect",
      values_to = "KAI"
    ) %>%
    rename(Fecha = `Fecha...1`, Long_date = `Fecha...2`) %>%
    mutate(Fecha = as.Date(Fecha, format = "%d/%m/%Y")) %>%
    filter(!is.na(KAI), !is.na(Fecha))

  census_date <- max(obs$Fecha)
  cat("Most recent census date:", format(census_date), "\n")

  # ---- Transects: sub-segment -> named-transect lookup ---------------------
  # Mirrors the transect layout used in function_rabbit_distances.R::get_sim_data
  transects <- vect(obs_transects)
  transects <- project(transects, crs(DT_template))

  area_trans <- data.frame(
    ID = 1:37,
    transect = c("Abalario", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Hinojos",
                 "Hinojos-Guadiamar", "Marismillas", "Muro", "Matochal", "Puntal",
                 "RBD-este", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea", "Sabinar-Mogea",
                 "Sabinar-Mogea", "RBD-este", "RBD-este", "Algaida-Sotos", "Algaida-Sotos",
                 "Algaida-Sotos", "Algaida-Sotos", "Coto del Rey", "Coto del Rey", "Coto del Rey",
                 "Coto del Rey", "Coto del Rey", "Coto del Rey", "Coto del Rey", "Hinojos",
                 "Hinojos", "Marismillas", "Marismillas", "Marismillas", "Puntal",
                 "Puntal", "Puntal")
  )
  area_trans$length <- (terra::perim(transects) / 1000)[area_trans$ID]

  extract_at_transects <- function(r_layer) {
    vals <- terra::extract(r_layer, transects, fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
    df <- area_trans
    df$value <- vals[df$ID]
    df %>%
      group_by(transect) %>%
      summarise(value = stats::weighted.mean(value, w = length, na.rm = TRUE), .groups = "drop")
  }

  nearest_layer <- function(r_stack, target_date) {
    t <- as.Date(time(r_stack))
    r_stack[[which.min(abs(t - target_date))]]
  }

  window_stat <- function(r_stack, target_date, months_back, fun) {
    t <- as.Date(time(r_stack))
    keep <- (t > (target_date %m-% months(months_back))) & (t <= target_date)
    if (!any(keep)) stop("No climate layers found in the ", months_back, "-month window before ", target_date)
    fun(r_stack[[keep]])
  }

  build_covariates_for_date <- function(d) {
    list(
      ndvi_mean = nearest_layer(r_ndvi_mean, d),
      ndvi_sd   = nearest_layer(r_ndvi_sd, d),
      flood     = nearest_layer(r_flood, d),
      height    = r_height,
      temp_6mo  = window_stat(r_temp, d, climate_window_months, function(x) mean(x, na.rm = TRUE)),
      precip_6mo= window_stat(r_precip, d, climate_window_months, function(x) sum(x, na.rm = TRUE))
    )
  }

  cat('\n-----------------------------------------------------------------')
  cat('\n         Building covariate table for all census dates')
  cat('\n-----------------------------------------------------------------\n')

  census_dates <- sort(unique(obs$Fecha))

  model_rows <- lapply(census_dates, function(d) {
    tryCatch({
      cov <- build_covariates_for_date(d)
      cov_list <- lapply(names(cov), function(nm) {
        extract_at_transects(cov[[nm]]) %>% rename(!!nm := value)
      })
      cov_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "transect"), cov_list)
      cov_df$Fecha <- d
      cov_df
    }, error = function(e) {
      cat("  Skipping", format(d), "- ", conditionMessage(e), "\n")
      NULL
    })
  }) %>% bind_rows()

  model_df <- obs %>%
    inner_join(model_rows, by = c("Fecha", "transect")) %>%
    filter(if_all(c(KAI, ndvi_mean, ndvi_sd, flood, height, temp_6mo, precip_6mo), ~ !is.na(.)))

  if (nrow(model_df) < 15) {
    stop("Not enough complete observation/covariate rows (", nrow(model_df), ") to fit a spatial model")
  }

  cat('\n-----------------------------------------------------------------')
  cat('\n         Fitting KAI ~ habitat + antecedent climate')
  cat('\n-----------------------------------------------------------------\n')

  # KAI (log1p-transformed to stay well-behaved with the zeros in the index)
  # regressed on habitat covariates and the mean/total climate of the months
  # leading up to each census.
  fit <- glm(log1p(KAI) ~ ndvi_mean + ndvi_sd + flood + height + temp_6mo + precip_6mo,
             data = model_df, family = gaussian())

  cat(sprintf("Fitted on %d observations (%d census dates)\n", nrow(model_df), length(unique(model_df$Fecha))))

  cat('\n-----------------------------------------------------------------')
  cat('\n         Predicting across space for the most recent census')
  cat('\n-----------------------------------------------------------------\n')

  pred_cov <- build_covariates_for_date(census_date)
  pred_stack <- rast(list(
    ndvi_mean  = pred_cov$ndvi_mean,
    ndvi_sd    = pred_cov$ndvi_sd,
    flood      = pred_cov$flood,
    height     = pred_cov$height,
    temp_6mo   = pred_cov$temp_6mo,
    precip_6mo = pred_cov$precip_6mo
  ))
  names(pred_stack) <- c("ndvi_mean", "ndvi_sd", "flood", "height", "temp_6mo", "precip_6mo")

  pred <- terra::predict(pred_stack, fit, se.fit = TRUE, index = 1:2)

  fit_link <- pred[[1]]
  se_link  <- pred[[2]]

  # Back-transform from log1p space; uncertainty via the delta method
  # (d/dx expm1(x) = exp(x)).
  mean_rast        <- expm1(fit_link)
  uncertainty_rast <- se_link * exp(fit_link)

  names(mean_rast)        <- "rabbit_population_mean"
  names(uncertainty_rast) <- "rabbit_population_uncertainty"

  mean_path        <- file.path(out_dir, mean_filename)
  uncertainty_path <- file.path(out_dir, uncertainty_filename)

  writeRaster(mean_rast, mean_path, overwrite = TRUE)
  writeRaster(uncertainty_rast, uncertainty_path, overwrite = TRUE)

  cat("Saved:", mean_path, "\n")
  cat("Saved:", uncertainty_path, "\n")

  invisible(list(
    census_date         = census_date,
    model               = fit,
    model_df            = model_df,
    mean_raster         = mean_rast,
    uncertainty_raster  = uncertainty_rast,
    mean_path           = mean_path,
    uncertainty_path    = uncertainty_path
  ))
}
