#
# Rabbit hindcast on the large (lynx) template: monthly rabbit density maps for the lynx prey maps
#
# Runs the calibrated rabbit model (posterior draws, as rabbit_backcast() does for the park) over
# the large template with the inputs made for it:
#   NDVI mean / within-cell SD  data/environmental_data/NDVI_large/harmonised/{mean,sd}.tif
#   flood (0/1)                 data/environmental_data/flood_monthly_500m_large.tif
#                               (environmental_data/create_large_flood_layer.R; includes the
#                               no-rabbit cells: towns, salines, water)
#   vegetation height           data/environmental_data/vegetation_height_500_large.tif
#                               (environmental_data/create_large_height_layer.R)
#   temperature, precipitation  the park files (ERA5-Land, <variable>_CDS), which cover the large area
#
# Outside the park the model is extrapolated: the K map parameters were calibrated on the park
# census only.
#
# Writes to out_dir:
#   rabbit_density_mean.tif     monthly mean density over the draws (time = first of the month)
#   rabbit_density_members.tif  every draw (layer names <YYYY_MM>_member_<i>)
#   sim_df.csv                  simulated KAI at the park census transects (check against census)
#   rhd_outbreaks.csv, rhd_monthly.csv (immunity mode)
#
# The first year is spin-up (the model starts at carrying capacity).
#
#   Rscript backcast/rabbit_hindcast_large.R [n_draws] [start_year] [end_year]

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(dplyr); library(tidyr)
  library(lubridate); library(sf); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
n_draws    <- if (length(args) >= 1) as.integer(args[1]) else 8
start_year <- if (length(args) >= 2) as.integer(args[2]) else 2005
end_year   <- if (length(args) >= 3) as.integer(args[3]) else 2025

ed <- file.path("data", "environmental_data")
template_file  <- file.path(ed, "template_raster_500_large.tif")
ndvi_mean_file <- file.path(ed, "NDVI_large", "harmonised", "mean.tif")
ndvi_sd_file   <- file.path(ed, "NDVI_large", "harmonised", "sd.tif")
flood_file     <- file.path(ed, "flood_monthly_500m_large.tif")
height_file    <- file.path(ed, "vegetation_height_500_large.tif")
temperature_dir <- file.path(ed, "temperature_CDS")
precip_file    <- file.path(ed, "precipitation_CDS")
census_dir     <- file.path("data", "rabbit", "Doñana_census_data_PacoCarro")
rates_file     <- file.path("data", "rabbit", "Rabbit_demography_values.txt")
model_exe      <- file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe")
# Latest immunity-mode calibration (park census 2017-2024)
posterior_file <- file.path("..", "Components", "output", "rabbit-calibration-immunity-pacc0.05", "posterior_samples.csv")
# RHD (classic) present from the start; RHD2 (GI.2) reached Doñana in 2011-2012. The run settings
# travel with them: a density map for every month (the default is March, June and September only).
rhd_settings_list <- function() modifyList(rhd_default_settings(),
                                           list(RHD2_arrival_year = 2011, RHD2_arrival_month = 10,
                                                create_all_month_maps = 1))
out_dir <- file.path("data", "rabbit_hindcast_large")

sapply(list.files(file.path("rabbit", "functions"), pattern = "\\.R$", full.names = TRUE), source)
source(file.path("backcast", "functions", "backcast_helpers.R"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir_clim <- file.path(out_dir, "clim_variables")
dir.create(out_dir_clim, showWarnings = FALSE)

template <- rast(template_file)
first_month <- as.Date(sprintf("%d-11-01", start_year - 1))
last_month  <- as.Date(sprintf("%d-12-01", end_year))
in_window <- function(r){ t <- as.Date(time(r)); r[[t >= first_month & t <= last_month]] }

# Climate. ERA5-Land is land-only, so coastal cells have no value: cells without a
# value take the nearest one (focal mean, repeated until full).
fill_edges <- function(r, mask){
  for (k in 1:10) {
    if (!any(is.na(values(mask(r[[1]], mask), mat = FALSE)) & !is.na(values(mask, mat = FALSE)))) break
    r <- cover(r, focal(r, 3, "mean", na.rm = TRUE))
  }
  mask(r, mask)
}
r_temp   <- in_window(load_monthly_stack(temperature_dir, template, pattern = "\\.nc$"))
r_precip <- in_window(load_monthly_stack(precip_file, template, pattern = "\\.nc$"))
r_temp   <- fill_edges(r_temp, template[[1]])
r_precip <- fill_edges(r_precip, template[[1]])

r_mean_ndvi <- in_window(read_monthly_ndvi(ndvi_mean_file))
r_sd_ndvi   <- in_window(read_monthly_ndvi(ndvi_sd_file))
r_flood     <- in_window(read_monthly_ndvi(flood_file))
r_height    <- rast(height_file)

checked <- check_time_span_list(list(temperature = r_temp, precipitation = r_precip,
                                     meanNDVI = r_mean_ndvi, sdNDVI = r_sd_ndvi),
                                start_year, end_year, fill_gaps = TRUE)

calculate_BM_cDM(temp = checked$temperature, K = TRUE, precip = checked$precipitation,
                 year_min = start_year, year_max = end_year,
                 result_dir = out_dir_clim, coord_rast = template)

Cov_list <- list(ndvi_mean = checked$meanNDVI, ndvi_sd = checked$sdNDVI,
                 flood = r_flood, height = r_height)

obs_df <- read_csv(file.path(census_dir, "KAI_Rabbit_Night_2024_v1.csv"), show_col_types = FALSE) %>%
  pivot_longer(cols = c("Coto del Rey", "Algaida-Sotos", "Sabinar-Mogea", "RBD-este", "Puntal",
                        "Marismillas", "Abalario", "Hinojos"),
               names_to = "transect", values_to = "KAI") %>%
  rename(Fecha = `Fecha...1`, Long_date = `Fecha...2`)

set.seed(20261108)
t0 <- Sys.time()
sim <- run_simulation(posterior_df = read.csv(posterior_file),
                      start_year = start_year, end_year = end_year,
                      N_samples = n_draws, compare_obs = TRUE,
                      Cov_list = Cov_list, DT_template = template,
                      obs_transects = normalizePath(file.path(census_dir, "Transect_oryctolagus.kml")),
                      obs_df = obs_df,
                      demograph_mod = normalizePath(model_exe),
                      rabbit_demographic_rates = normalizePath(rates_file),
                      demographic_climate = normalizePath(out_dir_clim),
                      rhd_settings = rhd_settings_list())
cat(sprintf("\n%d draws of %d-%d in %.1f min\n", n_draws, start_year, end_year,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Monthly maps in time order
ym <- regmatches(names(sim$pop_maps), regexpr("\\d{4}_\\d{1,2}", names(sim$pop_maps)))
dates <- as.Date(sprintf("%s-%02d-01", sub("_.*", "", ym), as.integer(sub(".*_", "", ym))))
o <- order(dates); maps <- sim$pop_maps[o]; dates <- dates[o]

mean_stack <- mask(rast(lapply(maps, `[[`, "mean")), template[[1]])
time(mean_stack) <- dates; names(mean_stack) <- format(dates, "%Y_%m")
writeRaster(mean_stack, file.path(out_dir, "rabbit_density_mean.tif"), overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE"))
members <- mask(rast(lapply(maps, `[[`, "members")), template[[1]])
names(members) <- paste(rep(format(dates, "%Y_%m"), each = n_draws), sprintf("member_%02d", seq_len(n_draws)), sep = "_")
time(members) <- rep(dates, each = n_draws)
writeRaster(members, file.path(out_dir, "rabbit_density_members.tif"), overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE"))
write_csv(sim$sim_df, file.path(out_dir, "sim_df.csv"))
if (!is.null(sim$rhd)) {
  write_csv(sim$rhd$outbreaks, file.path(out_dir, "rhd_outbreaks.csv"))
  write_csv(sim$rhd$monthly, file.path(out_dir, "rhd_monthly.csv"))
}

park <- rast(file.path(ed, "park_mask_500_large.tif"))
tot <- data.frame(date = dates,
                  park = global(mask(mean_stack, park), "sum", na.rm = TRUE)[, 1],
                  outside = global(mask(mean_stack, park, inverse = TRUE), "sum", na.rm = TRUE)[, 1])
write_csv(tot, file.path(out_dir, "total_abundance.csv"))
print(tot %>% mutate(year = year(date)) %>% group_by(year) %>%
        summarise(park = round(mean(park)), outside = round(mean(outside))), n = 30)
