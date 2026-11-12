library(terra)
library(lubridate)

# Small helpers for the backcast workflows. Need load_and_format_rasters() from
# rabbit/functions/miscellaneous.R.

# Monthly stack from several files or folders on the template grid (paths that don't exist are skipped). Unlike the calibration, the
# same month may appear in more than one file here (e.g. an older and a newer download of the
# same year), so duplicate months are averaged rather than summed.
load_monthly_stack <- function(paths, template, pattern = NULL, ...){
  paths <- paths[file.exists(paths)]          # e.g. backcast folders not made yet
  files <- unlist(lapply(paths, function(p){
    if(file_test("-d", p)) list.files(p, pattern = pattern, full.names = TRUE) else p
  }))
  load_and_format_rasters(files, template, temporal_agg = "month", temporal_fun = "mean", ...)
}

# Pre-aggregated monthly NDVI tif (e.g. NDVI_large/harmonised/mean.tif) with its time as Date.
# With a template (e.g. the park), the stack is cut to it: cropped and masked when the grids are
# aligned (the park grid is part of the large grid), otherwise resampled.
read_monthly_ndvi <- function(file, template = NULL){
  r <- rast(file)
  tr <- floor_date(as.Date(time(r)), "month")
  if(!is.null(template) && !compareGeom(r, template, stopOnError = FALSE)){
    r <- if(compareGeom(r, template, ext = FALSE, rowcol = FALSE, stopOnError = FALSE))
      crop(r, template) else resample(r, template)
    r <- mask(r, template[[1]])
  }
  time(r) <- tr
  r
}

# Months before the first layer of r (from start_date on) are filled with the average of the same
# calendar month over the first n_years years of data. Only meant for a spin-up period.
fill_leading_climatology <- function(r, start_date, n_years = 3){
  t <- as.Date(time(r))
  first <- min(t)
  start_date <- floor_date(as.Date(start_date), "month")
  if(first <= start_date) return(r)

  lead_months <- seq(start_date, first - months(1), by = "month")
  ref <- r[[t < first + years(n_years)]]
  ref_t <- as.Date(time(ref))

  filled <- lapply(lead_months, function(m){
    layer <- mean(ref[[month(ref_t) == month(m)]], na.rm = TRUE)
    time(layer) <- m
    layer
  })
  message(sprintf("Filled %d month(s) before %s with the %d-year monthly average",
                  length(lead_months), format(first, "%Y-%m"), n_years))

  out <- c(rast(filled), r)
  time(out) <- c(lead_months, t)
  out
}

