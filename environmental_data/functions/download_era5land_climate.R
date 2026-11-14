library(terra)
library(lubridate)
library(ecmwfr)

# Downloads ERA5-Land monthly means of 2 m temperature and total precipitation from the Copernicus
# Climate Data Store, on the native ERA5-Land grid (0.1 degree, ~9-11 km), without downscaling.
# The workflows bring them to their template when they load them (load_and_format_rasters in
# rabbit/functions/miscellaneous.R: project + bilinear resample).
#
# Output: one NetCDF per variable and calendar year,
#   <output_dir>/temperature_CDS/era5land_temperature_<yyyy>.nc     t2m, K
#   <output_dir>/precipitation_CDS/era5land_precipitation_<yyyy>.nc tp, m per day (monthly mean
#                                                                    of the daily total)
# with the time of each layer at the first of the month (00:00 UTC). The units are the ERA5-Land
# ones; the workflows convert them (K -> C, m/day -> mm).
#
# Years that already have all 12 months are skipped, so the function can be rerun to extend the
# series: the current year (and any year whose download was cut short) is downloaded again and
# its file replaced. ERA5-Land is published with ~2-3 months delay; the months not yet available
# are simply missing from the last file.
#
# Needs CDS credentials: API_Key (data/environmental_data/CDS_API_info.R) and the ERA5-Land
# licence accepted on the CDS website.
#
# area: c(north, west, south, east) in degrees. Use a margin of at least 0.2 degree around the
# template, so that the bilinear resampling has data at the template edges.
download_era5land_climate <- function(output_dir, area, api_key,
                                      start_year = 1984, end_year = year(Sys.Date()),
                                      variables = c(temperature   = "2m_temperature",
                                                    precipitation = "total_precipitation"),
                                      block_years = 10, retries = 3){

  wf_set_key(key = api_key)

  for(v in names(variables)){
    var_dir <- era5land_dir(output_dir, v)
    dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
    year_file <- function(y) file.path(var_dir, sprintf("era5land_%s_%d.nc", v, y))

    complete <- function(y){
      f <- year_file(y)
      file.exists(f) && nlyr(rast(f)) == 12
    }
    todo <- Filter(Negate(complete), start_year:end_year)
    if(length(todo) == 0){
      cat(v, ": all years", start_year, "-", end_year, "complete\n")
      next
    }

    # consecutive years in blocks, one CDS request per block
    blocks <- split(todo, cumsum(c(1, diff(todo) != 1)))
    blocks <- unlist(lapply(blocks, function(b) split(b, ceiling(seq_along(b) / block_years))),
                     recursive = FALSE)

    for(b in blocks){
      cat(v, ": downloading", min(b), "-", max(b), "\n")
      tmp_name <- sprintf("era5land_%s_%d_%d_download.nc", v, min(b), max(b))
      request <- list(
        dataset_short_name = "reanalysis-era5-land-monthly-means",
        product_type       = "monthly_averaged_reanalysis",
        variable           = variables[[v]],
        year               = as.character(b),
        month              = sprintf("%02d", 1:12),
        time               = "00:00",
        data_format        = "netcdf",
        download_format    = "unarchived",
        area               = area,
        target             = tmp_name
      )
      for(attempt in seq_len(retries)){
        ok <- tryCatch({ wf_request(request, path = tempdir(), time_out = 3600); TRUE },
                       error = function(e){ message("  attempt ", attempt, " failed: ", conditionMessage(e)); FALSE })
        if(ok) break
        if(attempt == retries) stop("CDS request for ", v, " ", min(b), "-", max(b), " failed ", retries, " times")
        Sys.sleep(60)
      }

      r <- rast(file.path(tempdir(), tmp_name))
      r <- r * 1                                   # into memory, so the download can be removed
      months <- floor_date(as.Date(time(r)), "month")
      varname <- if(v == "temperature") "t2m" else "tp"
      unit    <- if(v == "temperature") "K" else "m"

      for(y in b){
        i <- which(year(months) == y)
        if(length(i) == 0){
          warning(v, " ", y, ": no months available (yet)", call. = FALSE)
          next
        }
        ry <- r[[i]]
        time(ry) <- as.POSIXct(format(months[i]), tz = "UTC")
        names(ry) <- format(months[i], "%Y_%m")
        writeCDF(ry, year_file(y), varname = varname, unit = unit,
                 longname = if(v == "temperature") "2 metre temperature" else "total precipitation per day",
                 overwrite = TRUE)
        if(length(i) < 12) cat("  ", y, ": ", length(i), " months (to ", format(max(months[i]), "%Y-%m"), ")\n", sep = "")
      }
      rm(r); gc()
      unlink(file.path(tempdir(), tmp_name))
    }
  }
  invisible(NULL)
}

# Download area (c(north, west, south, east), degrees) for a template raster: its extent in
# WGS84 plus `margin` degrees, snapped outwards to the 0.1 degree ERA5-Land grid.
era5land_area <- function(template, margin = 0.3){
  e <- ext(project(as.polygons(ext(template), crs = crs(template)), "EPSG:4326"))
  c(ceiling((e$ymax + margin) * 10) / 10, floor((e$xmin - margin) * 10) / 10,
    floor((e$ymin - margin) * 10) / 10, ceiling((e$xmax + margin) * 10) / 10)
}

# Folder of one variable: <output_dir>/<variable>_CDS
era5land_dir <- function(output_dir, variable) file.path(output_dir, paste0(variable, "_CDS"))

# Monthly stack (one layer per month, sorted) of the per-year files of one variable
read_era5land_climate <- function(output_dir, variable = c("temperature", "precipitation")){
  variable <- match.arg(variable)
  var_dir <- era5land_dir(output_dir, variable)
  files <- list.files(var_dir, pattern = "^era5land_.*_\\d{4}\\.nc$", full.names = TRUE)
  if(length(files) == 0) stop("No ERA5-Land ", variable, " files in ", var_dir)
  r <- rast(lapply(files, rast))
  r[[order(time(r))]]
}
