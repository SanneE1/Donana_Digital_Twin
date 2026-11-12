#
# Rabbit backcast from 1985 on the large (lynx) template, with the rabbit components
#
# Prepares the inputs from the downloads, then runs the two LifeWatch components in
# ../Components:
#   donana-rabbit-calibration   ABC calibration on the park (cut out of the large inputs), simulated
#                               from 1985 so that the RHD arrival dates (1990, 2012) take effect,
#                               fitted to the park census 2005-2024
#   donana-rabbit-forecasting   the backcast: the posterior of the calibration run over the large
#                               template with the observed NDVI and climate (useNdviForecast and
#                               useClimateForecast FALSE), 1985 spin-up, 1986-01 to 2025-12
#
# Steps (each skips what is already on disk; force = TRUE redoes it):
#   climate      ERA5-Land monthly temperature and precipitation (native 0.1 degree grid, no
#                downscaling) for the large template, 1984 to 2025, in
#                data/environmental_data/<variable>_CDS. Needs CDS credentials
#                (data/environmental_data/CDS_API_info.R).
#   ndvi         monthly 500 m Landsat 5 and 7 products from the 1984+ downloads (the 2004+ products
#                are kept in <sensor>_2004), then Landsat adjusted to Sentinel-2 into
#                NDVI_large/harmonised_1984 (harmonised/ stays as it is)
#   flood        monthly LAST maps 1984+ from the WCS scenes (data/backcast/LAST_maps) on the large
#                template, then the large flood layer 1984-11 to 2025-12 (flood_monthly_500m_large_1984.tif)
#   inputs       one folder per component input, one file each: data/backcast/components/inputs_large
#                and inputs_park (the park cut out of the large inputs, same products)
#   calibration  runs donana-rabbit-calibration -> data/backcast/rabbit_backcast_1985/calibration
#   forecast     runs donana-rabbit-forecasting -> data/backcast/rabbit_backcast_1985/backcast,
#                plus the mean density over the runs, park/outside totals and the census agreement
#
# Run from the Donana_Digital_Twin project root:
#   Rscript backcast/Workflow_backcast_components.R                       # all steps
#   Rscript backcast/Workflow_backcast_components.R inputs calibration    # some steps
#
# With RABBIT_DOCKER=1 the components run in their Docker images (donana-rabbit-calibration:test,
# donana-rabbit-forecasting:test) with the Linux build of the model; otherwise with Rscript and the
# Windows build. Laptop memory: n_cluster 2 at most (see local machine notes).

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(dplyr); library(tidyr); library(lubridate); library(readr)
})

steps <- commandArgs(trailingOnly = TRUE)
all_steps <- c("climate", "ndvi", "flood", "inputs", "calibration", "forecast")
if (length(steps) == 0) steps <- all_steps
if (!all(steps %in% all_steps)) stop("Unknown step(s): ", paste(setdiff(steps, all_steps), collapse = ", "),
                                     ". Steps: ", paste(all_steps, collapse = ", "))
force <- FALSE

#-------------------------------------------------------------------------------
# Settings
#-------------------------------------------------------------------------------

first_month <- as.Date("1984-11-01")   # the climate of a model year starts in November before it
last_month  <- as.Date("2025-12-01")

# Calibration from 2003: two years before the census (2005-2024). RHD1 (1990) is then already
# present at the start (the model starts with the disease and part of the rabbits immune); RHD2
# (2012) arrives during the run.
calibration <- list(start_year = 2003, end_year = 2024,
                    nb_simul = 300, p_acc_min = 0.05, n_sim_samples = 20, n_cluster = 2,
                    rhdMode = "immunity")
backcast    <- list(startDate = "1986-01", spinupYears = 1, horizon = "",   # empty: to 2025-12
                    nRuns = 10, n_cluster = 2, uncertaintySources = "parameters,demographic",
                    rhdMode = "immunity")
# RHD arrivals, the same in both components (their defaults)
rhd_arrivals <- list(rhd1ArrivalYear = 1990, rhd1ArrivalMonth = 1, rhd2ArrivalYear = 2012, rhd2ArrivalMonth = 1)

ed <- file.path("data", "environmental_data")
large_template_file <- file.path(ed, "template_raster_500_large.tif")
park_template_file  <- file.path(ed, "template_raster_500.tif")
park_mask_file      <- file.path(ed, "park_mask_500_large.tif")

climate_dir   <- ed                                                   # <variable>_CDS/era5land_<variable>_<yyyy>.nc
ndvi_dir      <- file.path(ed, "NDVI_large")
ndvi_harm_dir <- file.path(ndvi_dir, "harmonised_1984")
last_scenes   <- file.path("data", "backcast", "LAST_maps")
last_monthly  <- file.path("data", "backcast", "LAST_monthly_floodmaps_large")
flood_file    <- file.path(ed, "flood_monthly_500m_large_1984.tif")
height_file   <- file.path(ed, "vegetation_height_500_large.tif")

census_dir    <- file.path("data", "rabbit", "Doñana_census_data_PacoCarro")
kai_file      <- file.path(census_dir, "KAI_Rabbit_Night_2024_v1.csv")
transect_file <- file.path(census_dir, "Transect_oryctolagus.kml")
rates_file    <- file.path("data", "rabbit", "Rabbit_demography_values.txt")

use_docker <- nzchar(Sys.getenv("RABBIT_DOCKER"))
model_exe  <- file.path("rabbit", "Pascal_program", "programs",
                        if (use_docker) "Rabbit_x86_64-linux" else "Rabbit_x86_64-win64.exe")
comp_root  <- normalizePath(file.path("..", "Components"), winslash = "/")
work_dir   <- file.path("data", "backcast", "components")
out_dir    <- file.path("data", "backcast", "rabbit_backcast_1985")

# BACKCAST_SMOKE=1: short test of the components (small calibration and backcast), written to
# rabbit_backcast_1985_smoke so that the real outputs are not touched
if (nzchar(Sys.getenv("BACKCAST_SMOKE"))) {
  calibration[c("nb_simul", "p_acc_min", "n_sim_samples")] <- list(68, 0.2, 5)
  backcast$nRuns <- 2
  work_dir <- file.path("data", "backcast", "components_smoke")
  out_dir  <- paste0(out_dir, "_smoke")
  cat("SMOKE TEST: outputs in", out_dir, "\n")
}
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sapply(list.files(file.path("rabbit", "functions"), pattern = "\\.R$", full.names = TRUE), source)
source(file.path("backcast", "functions", "backcast_helpers.R"))
source(file.path("environmental_data", "functions", "download_era5land_climate.R"))

template <- rast(large_template_file)
in_window <- function(r, from = first_month){ t <- as.Date(time(r)); r[[t >= from & t <= last_month]] }
months_covered <- function(r) format(range(as.Date(time(r))), "%Y-%m")
step_header <- function(s) cat("\n=================================================================",
                               "\n  ", s, " (", format(Sys.time(), "%H:%M"), ")",
                               "\n=================================================================\n", sep = "")

#-------------------------------------------------------------------------------
# 1. Climate (temperature and precipitation)
#-------------------------------------------------------------------------------

if ("climate" %in% steps) {
  step_header("Climate: ERA5-Land temperature and precipitation")
  source(file.path(ed, "CDS_API_info.R"))          # API_User, API_Key
  download_era5land_climate(output_dir = climate_dir, area = era5land_area(template), api_key = API_Key,
                            start_year = year(first_month), end_year = year(last_month))
  for (v in c("temperature", "precipitation"))
    cat(v, ":", months_covered(read_era5land_climate(climate_dir, v)), "\n")
}

#-------------------------------------------------------------------------------
# 2. NDVI: Landsat 5/7 products from 1984, adjusted to Sentinel-2
#-------------------------------------------------------------------------------

if ("ndvi" %in% steps && (force || !file.exists(file.path(ndvi_harm_dir, "mean.tif")))) {
  step_header("NDVI: monthly Landsat 5/7 products and harmonisation")
  source(file.path("environmental_data", "functions", "build_monthly_ndvi_500m.R"))
  for (s in c("landsat-5", "landsat-7")) {
    product_dir <- file.path(ndvi_dir, s)
    backup <- paste0(product_dir, "_2004")
    if (dir.exists(product_dir) && !dir.exists(backup)) {
      dir.create(backup)
      file.copy(list.files(product_dir, full.names = TRUE), backup)
      cat("Kept the 2004+ products of", s, "in", backup, "\n")
    }
    # raw monthly GEE downloads (ndvi_monthly_500m_gee.py); removed from disk on 2026-11-10 after the
    # products in NDVI_large were built, so they have to be downloaded again to rebuild
    raw_dir <- file.path(ed, "ndvi500_large", paste0("ndvi500_large_", s))
    if (!dir.exists(raw_dir)) stop("No raw NDVI downloads in ", raw_dir,
                                   ": download them again with backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py")
    build_monthly_ndvi_500m(raw_dir, template, product_dir)
  }
  e <- new.env()
  e$harmonise_area <- "large"
  e$harmonise_out_dir <- ndvi_harm_dir
  source(file.path("backcast", "harmonise_landsat_to_s2.R"), local = e)

  # Gaps the linear fill bridged: longest run of filled months per land cell, from 1985
  # (2026-11-09: median 3, 99% 4, max 7 months in 1985-2003; 2004-2015 median 5, max 7)
  filled <- in_window(read_monthly_ndvi(file.path(ndvi_harm_dir, "filled.tif")), as.Date("1985-01-01"))
  fv <- values(filled, mat = TRUE)[!is.na(values(template[[1]], mat = FALSE)), , drop = FALSE]
  run_max <- apply(fv, 1, function(x){ r <- rle(x %in% 1); max(c(0, r$lengths[r$values])) })
  q <- quantile(run_max, c(0.5, 0.9, 0.99, 1))
  cat("Longest gap filled per land cell (months) from 1985: median", q[1], "| 90%", q[2], "| 99%", q[3], "| max", q[4], "\n")
}

#-------------------------------------------------------------------------------
# 3. Flood layer from 1984
#-------------------------------------------------------------------------------

if ("flood" %in% steps && (force || !file.exists(flood_file))) {
  step_header("Flood: monthly LAST maps 1984+ and the large flood layer")
  source(file.path("backcast", "functions", "make_monthly_floodmaps.R"))
  scenes <- list.files(last_scenes, pattern = "\\.tif$", full.names = TRUE)
  bad <- scenes[file.size(scenes) < 2000]          # server error messages saved as .tif
  if (length(bad) > 0) cat("Skipped (not a map, < 2 KB):", basename(bad), "\n")
  make_monthly_floodmaps(setdiff(scenes, bad), template, last_monthly, min_year = 1984)

  e <- new.env()
  e$flood_months    <- seq(first_month, last_month, by = "month")
  # park months: the LAST maps the model was calibrated with (2000+) first, then the WCS ones
  e$flood_last_dirs <- c(file.path(ed, "LAST_monthly_floodmaps"), last_monthly)
  e$flood_out_file  <- flood_file
  source(file.path("environmental_data", "create_large_flood_layer.R"), local = e)
}

#-------------------------------------------------------------------------------
# 4. Component inputs: one folder per input, one file per folder
#-------------------------------------------------------------------------------

# The components load every file in an input folder, so each gets its own folder with one GeoTIFF
# (time on the layers; no .aux.xml next to it)
write_input <- function(r, dir, name, timed = TRUE){
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  f <- file.path(dir, name)
  writeRaster(r, f, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
  unlink(list.files(dir, pattern = "\\.aux\\.xml$", full.names = TRUE))
  if (timed && all(is.na(time(rast(f))))) stop("No time information in ", f)
  cat(sprintf("  %-55s %3d layer(s)%s\n", f, nlyr(r), if (timed) paste(",", paste(months_covered(r), collapse = " to ")) else ""))
}

# Nearest-value fill for cells at the edge without a climate value
fill_edges <- function(r, mask){
  for (k in 1:10) {
    if (!any(is.na(values(mask(r[[1]], mask), mat = FALSE)) & !is.na(values(mask, mat = FALSE)))) break
    r <- cover(r, focal(r, 3, "mean", na.rm = TRUE))
  }
  mask(r, mask)
}

inputs_large <- file.path("data", "backcast", "components", "inputs_large")   # also for the smoke test
inputs_park  <- file.path("data", "backcast", "components", "inputs_park")

if ("inputs" %in% steps && (force || !dir.exists(inputs_park))) {
  step_header("Component inputs")
  vars <- list(
    temperature   = fill_edges(in_window(load_monthly_stack(era5land_dir(climate_dir, "temperature"), template,
                                                            pattern = "\\.nc$")), template[[1]]),
    precipitation = fill_edges(in_window(load_monthly_stack(era5land_dir(climate_dir, "precipitation"), template,
                                                            pattern = "\\.nc$")), template[[1]]),
    ndvi_mean     = in_window(read_monthly_ndvi(file.path(ndvi_harm_dir, "mean.tif"))),
    ndvi_sd       = in_window(read_monthly_ndvi(file.path(ndvi_harm_dir, "sd.tif"))),
    flood         = in_window(read_monthly_ndvi(flood_file)))

  # every month present?
  needed <- format(seq(first_month, last_month, by = "month"), "%Y-%m")
  for (v in names(vars)) {
    from <- if (v %in% c("temperature", "precipitation")) 1 else 3     # NDVI and flood from January
    miss <- setdiff(needed[from:length(needed)], format(as.Date(time(vars[[v]])), "%Y-%m"))
    if (length(miss) > 0) stop(v, " is missing ", length(miss), " month(s): ", paste(head(miss, 12), collapse = ", "))
    cells_na <- global(is.na(vars[[v]][[1]]) & !is.na(template[[1]]), "sum", na.rm = TRUE)[1, 1]
    if (cells_na > 0) warning(v, ": ", cells_na, " template cells without a value in the first month", call. = FALSE)
  }
  height <- rast(height_file)

  park <- rast(park_template_file)
  to_park <- function(r) mask(crop(r, park), park[[1]])

  for (v in names(vars)) {
    write_input(vars[[v]], file.path(inputs_large, v), paste0(v, ".tif"))
    write_input(to_park(vars[[v]]), file.path(inputs_park, v), paste0(v, ".tif"))
  }
  write_input(height, file.path(inputs_large, "height"), "vegetation_height.tif", timed = FALSE)
  write_input(to_park(height), file.path(inputs_park, "height"), "vegetation_height.tif", timed = FALSE)
}

#-------------------------------------------------------------------------------
# Running a component
#-------------------------------------------------------------------------------

# inputs: named list, input name -> list(from = file or folder, to = path under /mnt/inputs/<name>/)
stage_component <- function(sandbox, inputs){
  unlink(sandbox, recursive = TRUE)
  for (nm in names(inputs)) {
    to <- file.path(sandbox, "mnt", "inputs", nm, inputs[[nm]]$to)
    from <- inputs[[nm]]$from
    if (dir.exists(from)) {
      dir.create(to, recursive = TRUE, showWarnings = FALSE)
      stopifnot(all(file.copy(list.files(from, full.names = TRUE), to)))
    } else {
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      stopifnot(file.exists(from), file.copy(from, to, overwrite = TRUE))
    }
  }
  dir.create(file.path(sandbox, "mnt", "outputs"), recursive = TRUE, showWarnings = FALSE)
  normalizePath(sandbox, winslash = "/")
}

run_component <- function(component, sandbox, params){
  args <- paste0("--", names(params), "=", unlist(params))
  cat("Running", component, "with", paste(args, collapse = " "), "\n")
  t0 <- Sys.time()
  status <- if (use_docker) {
    system2("docker", c("run", "--rm",
                        "-v", paste0(file.path(sandbox, "mnt", "inputs"), ":/mnt/inputs"),
                        "-v", paste0(file.path(sandbox, "mnt", "outputs"), ":/mnt/outputs"),
                        paste0(component, ":test"), args))
  } else {
    old <- setwd(file.path(comp_root, component))
    on.exit(setwd(old), add = TRUE)
    Sys.setenv(LW_SIM_ROOT = sandbox)
    system2("Rscript", c("--vanilla", "main.r", args))
  }
  cat(sprintf("%s finished in %.0f min (status %s)\n", component,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), status))
  if (status != 0) stop(component, " failed (status ", status, "); the sandbox is ", sandbox)
  file.path(sandbox, "mnt", "outputs")
}

copy_outputs <- function(from, to){
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(from, full.names = TRUE), to, recursive = TRUE, overwrite = TRUE)
  cat("Outputs in", to, "\n")
}

common_inputs <- function(inputs_dir, template_file) list(
  templateRaster   = list(from = template_file, to = "template_raster_500.tif"),
  obsTransects     = list(from = transect_file, to = "Transect_oryctolagus.kml"),
  obsKAI           = list(from = kai_file, to = "KAI_Rabbit_Night.csv"),
  demographicRates = list(from = rates_file, to = "Rabbit_demography_values.txt"),
  # main.r gives a Windows build its .exe extension
  demographicModel = list(from = model_exe, to = "Rabbit_x86_64-linux"),
  HeightFile       = list(from = file.path(inputs_dir, "height", "vegetation_height.tif"), to = "vegetation_height.tif"),
  TempDir          = list(from = file.path(inputs_dir, "temperature"), to = ""),
  NdviMaxDir       = list(from = file.path(inputs_dir, "ndvi_mean"), to = ""),
  NdviSdDir        = list(from = file.path(inputs_dir, "ndvi_sd"), to = ""),
  FloodDir         = list(from = file.path(inputs_dir, "flood"), to = ""))

posterior_file <- file.path(out_dir, "calibration", "posterior_samples.csv")

#-------------------------------------------------------------------------------
# 5. Calibration (park, 1985-2024)
#-------------------------------------------------------------------------------

if ("calibration" %in% steps && (force || !file.exists(posterior_file))) {
  step_header("Calibration: donana-rabbit-calibration on the park, 1985-2024")
  inputs <- c(common_inputs(inputs_park, park_template_file),
              list(PrecipFile = list(from = file.path(inputs_park, "precipitation"), to = "Precipitation")))
  sandbox <- stage_component(file.path(work_dir, "calibration"), inputs)
  outputs <- run_component("donana-rabbit-calibration", sandbox, c(calibration, rhd_arrivals))
  copy_outputs(outputs, file.path(out_dir, "calibration"))
}

#-------------------------------------------------------------------------------
# 6. Backcast (large template, observed NDVI and climate)
#-------------------------------------------------------------------------------

if ("forecast" %in% steps) {
  step_header("Backcast: donana-rabbit-forecasting on the large template")
  if (!file.exists(posterior_file)) stop("No posterior yet: run the calibration step first (", posterior_file, ")")
  inputs <- c(common_inputs(inputs_large, large_template_file),
              list(PrecipDir        = list(from = file.path(inputs_large, "precipitation"), to = ""),
                   posteriorSamples = list(from = posterior_file, to = "posterior_samples.csv")))
  sandbox <- stage_component(file.path(work_dir, "forecasting"), inputs)
  params <- c(backcast, rhd_arrivals,
              list(useNdviForecast = "FALSE", useClimateForecast = "FALSE", runMaps = "TRUE"))
  params <- params[vapply(params, function(p) nzchar(as.character(p)), logical(1))]
  outputs <- run_component("donana-rabbit-forecasting", sandbox, params)
  bc_dir <- file.path(out_dir, "backcast")
  copy_outputs(outputs, bc_dir)

  # Mean monthly density over the runs, and totals in and outside the park
  run_files <- list.files(file.path(bc_dir, "run_density_maps"), pattern = "\\.tif$", full.names = TRUE)
  runs <- lapply(run_files, rast)
  mean_stack <- Reduce(`+`, runs) / length(runs)
  time(mean_stack) <- time(runs[[1]]); names(mean_stack) <- names(runs[[1]])
  writeRaster(mean_stack, file.path(out_dir, "rabbit_density_mean.tif"), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

  park <- rast(park_mask_file)
  tot <- data.frame(date    = as.Date(time(mean_stack)),
                    park    = global(mask(mean_stack, park, maskvalues = c(0, NA)), "sum", na.rm = TRUE)[, 1],
                    outside = global(mask(mean_stack, park, maskvalues = 1), "sum", na.rm = TRUE)[, 1])
  write_csv(tot, file.path(out_dir, "total_abundance.csv"))
  print(tot %>% mutate(year = year(date)) %>% group_by(year) %>%
          summarise(park = round(mean(park)), outside = round(mean(outside))), n = 50)

  # Agreement with the park census: ensemble-mean KAI against the observed KAI per transect
  obs <- read_csv(kai_file, show_col_types = FALSE, name_repair = "unique_quiet") %>%
    pivot_longer(c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea", "RBD-este", "Puntal",
                   "Marismillas", "Abalario", "Hinojos"), names_to = "transect", values_to = "KAI_obs") %>%
    mutate(Fecha = as.Date(`Fecha...1`, format = "%d/%m/%Y")) %>% select(Fecha, transect, KAI_obs)
  sim <- read_csv(file.path(bc_dir, "forecast_kai_summary.csv"), show_col_types = FALSE) %>%
    filter(scenario == "all") %>% select(Fecha, transect, KAI_sim = mean)
  agree <- inner_join(sim, obs, by = c("Fecha", "transect")) %>% filter(!is.na(KAI_obs)) %>%
    group_by(transect) %>%
    summarise(n = n(), spearman = suppressWarnings(cor(KAI_sim, KAI_obs, method = "spearman")), .groups = "drop")
  write_csv(agree, file.path(out_dir, "agreement_with_census.csv"))
  print(agree)
}
