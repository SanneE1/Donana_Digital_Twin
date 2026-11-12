#
# Download monthly ERA5-Land temperature and precipitation, 1984 to the latest month available
#
# Area: the large (lynx) template plus a 0.3 degree margin, which includes the park template.
# Output (one file per variable and year, native 0.1 degree grid, no downscaling):
#   data/environmental_data/temperature_CDS/era5land_temperature_<yyyy>.nc     (K)
#   data/environmental_data/precipitation_CDS/era5land_precipitation_<yyyy>.nc (m/day)
# Rerun to extend the series: complete years are skipped, the current year is refreshed.
#
# Needs data/environmental_data/CDS_API_info.R (API_User, API_Key) and the ecmwfr package.
#
# Run from the Donana_Digital_Twin project root:
#   Rscript environmental_data/download_climate.R            # 1984 to now
#   Rscript environmental_data/download_climate.R 2020 2021  # only these years

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))

data_dir    <- file.path("data", "environmental_data")
climate_dir <- data_dir                          # files go to <climate_dir>/<variable>_CDS/
template    <- terra::rast(file.path(data_dir, "template_raster_500_large.tif"))

source(file.path(data_dir, "CDS_API_info.R"))   # API_User, API_Key
source(file.path("environmental_data", "functions", "download_era5land_climate.R"))

args <- as.integer(commandArgs(trailingOnly = TRUE))
start_year <- if(length(args) >= 1) args[1] else 1984
end_year   <- if(length(args) >= 2) args[2] else lubridate::year(Sys.Date())

area <- era5land_area(template)
cat("Area (N, W, S, E):", area, "\n")

download_era5land_climate(output_dir = climate_dir, area = area, api_key = API_Key,
                          start_year = start_year, end_year = end_year)

# Coverage
for(v in c("temperature", "precipitation")){
  r <- read_era5land_climate(climate_dir, v)
  m <- as.Date(time(r))
  expected <- seq(min(m), max(m), by = "month")
  cat(sprintf("%-13s %s to %s, %d months, %d missing\n", v, format(min(m), "%Y-%m"),
              format(max(m), "%Y-%m"), length(m), sum(!expected %in% m)))
}
