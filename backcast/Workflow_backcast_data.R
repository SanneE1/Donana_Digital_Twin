#
# Backcast: download and prepare environmental data from the earliest possible date up to now
#
# Builds on the forecast workflow's corrected inputs, which are only read here:
#   - NDVI: data/environmental_data/NDVI_large/harmonised (500 m monthly mean/SD, Landsat adjusted
#     to Sentinel-2, gap-filled; the park is cut out of it). The park Sentinel-2 product
#     monthly_ndvi_v2 and the extent shapefile Donana_DT_template_extent were removed (2026-11-10).
# downloads the ERA5-Land temperature and precipitation from the climate start on (shared with the
# forecast workflow, data/environmental_data/<variable>_CDS), and adds on top, in data/backcast:
#   - monthly flood maps from 1984 (the forecast workflow starts in 2000)
#   - Landsat 5/7/8 NDVI per sensor (phase 2), adjusted afterwards by harmonise_landsat_to_s2.R
#
# phase = 1: backcast from 2015 (Sentinel-2). Nov/Dec 2014 climate is needed for the rabbit
#            climate variables, and is in the forecast workflow's files already.
# phase = 2: backcast from 1985 (Landsat before July 2015).
#
# Run from the Donana_Digital_Twin project root. Every step skips what is already on disk.

library(terra)
library(lubridate)
library(ncdf4)

phase <- 1

# Set paths---------------------------------------------------------------------
data_dir      = file.path("data", "environmental_data")   # forecast workflow data, read only
bc_dir        = file.path("data", "backcast")
aoi_path      = file.path(data_dir, "Donana_DT_template_extent", "DT_template_extent.shp")
template_path = file.path(data_dir, "template_raster_500.tif")
dem_path      = file.path(data_dir, "Copernicus_GLO90_Europe_250m.tif")

backcast_start <- if(phase == 1) as.Date("2015-01-01") else as.Date("1985-01-01")
climate_start  <- backcast_start - months(2)

dir.create(bc_dir, recursive = T, showWarnings = F)

# load functions----------------------------------------------------------------
# this file just assigns two objects: API_User and API_Key but stored in data to keep it private
if(file.exists(file.path(data_dir, "CDS_API_info.R"))){
  source(file.path(data_dir, "CDS_API_info.R"))
} else {
  API_User = NULL
  API_Key = NULL
}

source(file.path("rabbit", "functions", "miscellaneous.R"))
source(file.path("environmental_data", "functions", "download_era5land_climate.R"))
source(file.path("environmental_data", "functions", "build_monthly_ndvi.R"))
sapply(list.files(file.path("backcast", "functions"), pattern = "\\.R$", full.names = T), source)

DT_template <- rast(template_path)

# The extent shapefile is only needed for the Landsat downloads of phase 2. It was made by
# environmental_data/rebuild_monthly_ndvi.R (retired) and removed on 2026-11-10.
if(phase == 2 && !file.exists(aoi_path)){
  stop(aoi_path, " is needed for the Landsat downloads of phase 2 but was removed. The large-area ",
       "Landsat NDVI now comes from Google Earth Engine (backcast/ndvi_pipeline/ndvi_monthly_500m_gee.py)")
}

# Same clean python environment as the forecast workflow (see there)
env_dir <- "python_env"
python_exe <- file.path(env_dir, "Scripts", "python.exe")

if (!file.exists(python_exe)) {
  system(sprintf('python -m venv "%s"', env_dir))
  system(sprintf('"%s" -m pip install --upgrade pip', python_exe))
  system2(python_exe, args = c("-m", "pip", "install", "-r", "backcast/ndvi_pipeline/requirements.txt"))
  system2(python_exe, args = c("-m", "pip", "install", "owslib"))
}


#-------------------------------------------------------------------------------
# Temperature and precipitation: ERA5-Land monthly means, 1984 to now
#-------------------------------------------------------------------------------

# Shared with the forecast workflow: data/environmental_data/<variable>_CDS, one file per
# year on the native 0.1 degree grid; only missing or incomplete years are downloaded.
climate_dir <- data_dir
download_era5land_climate(output_dir = climate_dir,
                          area       = era5land_area(rast(file.path(data_dir, "template_raster_500_large.tif"))),
                          api_key    = API_Key,
                          start_year = year(climate_start))


#-------------------------------------------------------------------------------
# Flooding: LAST flood maps from the ICTS WCS (from 1984) + manual downloads of the recent years
#-------------------------------------------------------------------------------

system2(python_exe,
        args = c("backcast/functions/download_flooding_wcs_backcast.py", file.path(bc_dir, "LAST_maps")))

fld_files <- c(list.files(file.path(data_dir, "Inundacion"), full.names = T),
               list.files(file.path(bc_dir, "LAST_maps"), full.names = T))

make_monthly_floodmaps(fld_files = fld_files,
                       template  = DT_template,
                       out_dir   = file.path(bc_dir, "monthly_floodmaps"),
                       min_year  = year(climate_start))


#-------------------------------------------------------------------------------
# NDVI - Landsat (phase 2)
#-------------------------------------------------------------------------------

# One series per sensor, so each can be adjusted to Sentinel-2 separately:
#   Landsat 5: the only sensor before 1999
#   Landsat 7: overlaps Landsat 5 (1999-2011) and Sentinel-2 (2015-2022)
#   Landsat 8: overlaps Sentinel-2 (2015 onwards)
# Built with the same monthly processing as Sentinel-2 (build_monthly_ndvi); only the raw
# (unfilled) products are used by harmonise_landsat_to_s2.R.
landsat_years <- list("landsat-5" = 1984:2011,
                      "landsat-7" = 1999:2022,
                      "landsat-8" = 2013:2025)

if(phase == 2){
  ls_dir <- file.path(bc_dir, "NDVI_landsat")
  dir.create(ls_dir, recursive = T, showWarnings = F)

  for(platform in names(landsat_years)){
    for(yr in landsat_years[[platform]]){
      out_name <- sprintf("ndvi_%s_%d", platform, yr)
      if(file.exists(file.path(ls_dir, paste0(out_name, ".nc")))) next
      status <- system2(python_exe,
                        args = c("backcast/ndvi_pipeline/pipeline.py",
                                 "--source", "landsat", "--platforms", platform,
                                 "--dates", sprintf("%d-01-01/%d-01-01", yr, yr + 1),
                                 "--output-dir", ls_dir, "--output-name", out_name,
                                 "--shape_file", aoi_path))
      if(status != 0) warning("Landsat download failed for ", out_name, call. = F)
    }

    build_monthly_ndvi(list.files(ls_dir, pattern = paste0("^ndvi_", platform, "_\\d{4}\\.nc$"), full.names = T),
                       DT_template, file.path(bc_dir, paste0("monthly_ndvi_", platform)))
  }
}


#-------------------------------------------------------------------------------
# Coverage report
#-------------------------------------------------------------------------------

report_coverage <- function(name, r){
  t <- floor_date(as.Date(time(r)), "month")
  expected <- seq(floor_date(min(t), "month"), floor_date(max(t), "month"), by = "month")
  cat(sprintf("%-14s %s to %s, %d months, %d missing, %d duplicated",
              name, format(min(t), "%Y-%m"), format(max(t), "%Y-%m"), length(unique(t)),
              sum(!expected %in% t), sum(duplicated(t))), fill = T)
}

report_coverage("NDVI", rast(file.path(data_dir, "NDVI_large", "harmonised", "mean.tif")))
report_coverage("Flooding", load_and_format_rasters(file.path(bc_dir, "monthly_floodmaps"), DT_template))
report_coverage("Temperature", read_era5land_climate(climate_dir, "temperature"))
report_coverage("Precipitation", read_era5land_climate(climate_dir, "precipitation"))
