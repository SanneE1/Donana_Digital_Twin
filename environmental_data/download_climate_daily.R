#
# Download DAILY ERA5-Land temperature and precipitation, 1985 to the latest day available.
# The observed climate for the climate delta correction component (bias correction of the
# AdapteCCa projections, Components/climate-delta-correction).
#
# Area: the large (lynx) template plus a 0.3 degree margin, as the monthly download.
# Output (one file per variable and year, native 0.1 degree grid):
#   data/environmental_data/temperature_daily_CDS/era5land_temperature_daily_<yyyy>.nc     (K)
#   data/environmental_data/precipitation_daily_CDS/era5land_precipitation_daily_<yyyy>.nc (m/day)
# Rerun to extend the series: complete years are skipped. ~1-2 min per variable and year.
#
# Needs data/environmental_data/CDS_API_info.R (API_User, API_Key) and the ecmwfr package.
#
# Run from the Donana_Digital_Twin project root:
#   Rscript environmental_data/download_climate_daily.R            # 1985 to now
#   Rscript environmental_data/download_climate_daily.R 2000 2001  # only these years

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))

data_dir <- file.path("data", "environmental_data")
template <- terra::rast(file.path(data_dir, "template_raster_500_large.tif"))

source(file.path(data_dir, "CDS_API_info.R"))   # API_User, API_Key
source(file.path("environmental_data", "functions", "download_era5land_climate.R"))   # era5land_area()
source(file.path("environmental_data", "functions", "download_era5land_daily.R"))

args <- as.integer(commandArgs(trailingOnly = TRUE))
start_year <- if(length(args) >= 1) args[1] else 1985
end_year   <- if(length(args) >= 2) args[2] else lubridate::year(Sys.Date())

area <- era5land_area(template)
cat("Area (N, W, S, E):", area, "\n")

download_era5land_daily(output_dir = data_dir, area = area, api_key = API_Key,
                        start_year = start_year, end_year = end_year)

# Coverage
for(v in c("temperature", "precipitation")){
  files <- list.files(era5land_daily_dir(data_dir, v), pattern = "\\.nc$", full.names = TRUE)
  if(length(files) == 0) next
  d <- as.Date(unlist(lapply(files, function(f) as.character(as.Date(terra::time(terra::rast(f)))))))
  expected <- seq(min(d), max(d), by = "day")
  cat(sprintf("%-13s %s to %s, %d days, %d missing\n", v, min(d), max(d), length(d),
              sum(!expected %in% d)))
}
