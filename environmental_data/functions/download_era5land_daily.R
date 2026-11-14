library(terra)
library(lubridate)
library(ecmwfr)

# Downloads DAILY ERA5-Land 2 m temperature and total precipitation from the Copernicus Climate
# Data Store, on the native ERA5-Land grid (0.1 degree), e.g. as the observed climate of the
# climate delta correction component (Components/climate-delta-correction).
#
#   temperature:   daily mean of the 24 hourly values of each UTC day (dataset
#                  reanalysis-era5-land, averaged here), K
#   precipitation: daily total, m per day. ERA5-Land precipitation is accumulated from 00 UTC,
#                  so the 00 UTC value of day D+1 (dataset reanalysis-era5-land) is the total of
#                  day D; the dates are shifted back by one day accordingly.
#
# Output: one NetCDF per variable and calendar year,
#   <output_dir>/temperature_daily_CDS/era5land_temperature_daily_<yyyy>.nc     t2m, K
#   <output_dir>/precipitation_daily_CDS/era5land_precipitation_daily_<yyyy>.nc tp,  m/day
# Years that already have all their days are skipped, so the function can be rerun to extend
# the series (ERA5-Land is published with a few days' delay; a year still running is replaced).
#
# Needs CDS credentials (API_Key in data/environmental_data/CDS_API_info.R) and the ERA5-Land
# licence accepted on the CDS website.
#
# area: c(north, west, south, east) in degrees; era5land_area() in download_era5land_climate.R
# gives it for a template raster.
download_era5land_daily <- function(output_dir, area, api_key,
                                    start_year = 1985, end_year = year(Sys.Date()),
                                    variables = c("temperature", "precipitation"),
                                    retries = 5){

  wf_set_key(key = api_key)
  variables <- match.arg(variables, several.ok = TRUE)

  for(v in variables){
    var_dir <- era5land_daily_dir(output_dir, v)
    dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
    year_file <- function(y) file.path(var_dir, sprintf("era5land_%s_daily_%d.nc", v, y))
    n_days <- function(y) as.integer(as.Date(sprintf("%d-12-31", y)) - as.Date(sprintf("%d-01-01", y))) + 1L

    complete <- function(y){
      f <- year_file(y)
      file.exists(f) && nlyr(rast(f)) == n_days(y)
    }
    todo <- Filter(Negate(complete), start_year:end_year)
    if(length(todo) == 0){
      cat(v, ": all years", start_year, "-", end_year, "complete\n")
      next
    }

    for(y in todo){
      cat(v, ": downloading", y, "\n")
      r <- if(v == "temperature") request_t2m_daily(y, area, retries) else request_tp_daily(y, area, retries)
      if(is.null(r)){
        warning(v, " ", y, ": no days available (yet)", call. = FALSE)
        next
      }
      d <- as.Date(time(r))
      time(r) <- d
      names(r) <- format(d, "%Y-%m-%d")
      writeCDF(r, year_file(y), varname = if(v == "temperature") "t2m" else "tp",
               unit = if(v == "temperature") "K" else "m",
               longname = if(v == "temperature") "2 metre temperature, daily mean" else "total precipitation per day",
               compression = 4, overwrite = TRUE)
      if(nlyr(r) < n_days(y)) cat("  ", y, ": ", nlyr(r), " days (to ", format(max(d)), ")\n", sep = "")
    }
  }
  invisible(NULL)
}

# Folder of one variable: <output_dir>/<variable>_daily_CDS
era5land_daily_dir <- function(output_dir, variable) file.path(output_dir, paste0(variable, "_daily_CDS"))

# One CDS request with retries; returns the downloaded file, or NULL when the CDS has no data.
# The CDS sometimes fails every request for a while (3 failures in 13 min, November 2026), so
# the attempts are 5 minutes apart.
cds_request <- function(request, retries){
  for(attempt in seq_len(retries)){
    f <- tryCatch(wf_request(request, path = tempdir(), time_out = 3 * 3600),
                  error = function(e){ message("  attempt ", attempt, " failed: ", conditionMessage(e)); NULL })
    if(!is.null(f) && file.exists(f)) return(f)
    if(attempt < retries) Sys.sleep(300)
  }
  stop("CDS request ", request$target, " failed ", retries, " times")
}

# months of year y that have (started to) happen
months_of <- function(y){
  last <- if(y < year(Sys.Date())) 12 else month(Sys.Date())
  sprintf("%02d", seq_len(last))
}

# Daily mean temperature of year y: the hourly values (00-23 UTC) from reanalysis-era5-land,
# one request per month, averaged per UTC day here. (The CDS dataset
# derived-era5-land-daily-statistics gives the same daily means, but its requests sat in the
# CDS queue for hours, November 2026.) Days with fewer than 24 hours (the last days
# published) are dropped.
request_t2m_daily <- function(y, area, retries){
  rs <- lapply(months_of(y), function(m){
    cat("   month", m, "\n")
    req <- list(dataset_short_name = "reanalysis-era5-land", variable = "2m_temperature",
                year = as.character(y), month = m, day = sprintf("%02d", 1:31),
                time = sprintf("%02d:00", 0:23),
                data_format = "netcdf", download_format = "unarchived", area = area,
                target = sprintf("era5land_t2m_hourly_%d_%s_download.nc", y, m))
    f <- tryCatch(cds_request(req, retries), error = function(e){
      # the current month may not be available yet
      if(y == year(Sys.Date()) && as.integer(m) == month(Sys.Date())) NULL else stop(e)
    })
    if(is.null(f)) return(NULL)
    h <- rast(f) * 1
    unlink(f)
    day <- format(time(h), "%Y-%m-%d", tz = "UTC")
    full <- names(which(table(day) == 24))
    h <- h[[which(day %in% full)]]
    if(nlyr(h) == 0) return(NULL)
    r <- tapp(h, index = day[day %in% full], fun = mean)
    time(r) <- as.Date(sort(full))
    r
  })
  rs <- Filter(Negate(is.null), rs)
  if(length(rs) == 0) return(NULL)
  rast(rs)
}

# Daily precipitation totals of year y: the 00 UTC accumulations from 2 January of y to
# 1 January of y + 1 (two requests), dated one day earlier
request_tp_daily <- function(y, area, retries){
  get <- function(yy, months, days, tag){
    req <- list(dataset_short_name = "reanalysis-era5-land", variable = "total_precipitation",
                year = as.character(yy), month = months, day = days, time = "00:00",
                data_format = "netcdf", download_format = "unarchived", area = area,
                target = sprintf("era5land_tp_00utc_%s_download.nc", tag))
    f <- cds_request(req, retries)
    r <- rast(f) * 1
    unlink(f)
    r
  }
  r <- get(y, months_of(y), sprintf("%02d", 1:31), y)
  if(as.Date(sprintf("%d-01-01", y + 1)) <= Sys.Date()){
    nxt <- tryCatch(get(y + 1, "01", "01", paste0(y + 1, "_0101")), error = function(e) NULL)
    if(!is.null(nxt)) r <- c(r, nxt)
  }
  d <- as.Date(time(r)) - 1
  keep <- which(year(d) == y)
  if(length(keep) == 0) return(NULL)
  r <- r[[keep]]
  time(r) <- d[keep]
  r
}
