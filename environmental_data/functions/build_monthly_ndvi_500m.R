library(terra)
library(lubridate)

# Monthly NDVI model inputs from the per-month 500 m files of backcast/ndvi_pipeline/ndvi_monthly_500m.py
# or ndvi_monthly_500m_gee.py (Google Earth Engine)
#
# Same products and file names as build_monthly_ndvi() (which works from the 10 m files of
# pipeline.py), so either can be used as model input:
#   raw_mean.tif, raw_sd.tif, valid_fraction.tif   per month, before masking and filling
#   mean.tif, sd.tif                               cells with valid fraction < min_valid_frac masked,
#                                                  then gap-filled over time: the model input
#   filled.tif                                     1 where mean/sd were filled, 0 where observed
#
# ndvi500_dir: folder with ndvi500_<sensor>_<YYYY-MM>.tif files (months without scenes have a
# .none file and become gaps). fill_ndvi_gaps() is the one of build_monthly_ndvi.R.

source(file.path("environmental_data", "functions", "build_monthly_ndvi.R"))

# Puts a file of ndvi_monthly_500m.py (on the template grid) or of ndvi_monthly_500m_gee.py (template
# grid shifted by < 10 m, WGS 84 / UTM 29N, -9999 = no data) on the template grid
on_template <- function(r, template){
  r <- classify(r, cbind(-9999, NA))
  if(compareGeom(r, template, stopOnError = FALSE)) return(r)
  shift <- c(xmin(template) - xmin(r), ymax(template) - ymax(r))
  if(any(abs(shift) >= 10) || abs(res(r)[1] - res(template)[1]) > 1e-6 ||
     !grepl("UTM zone 29N", crs(r, describe = TRUE)$name)){
    stop("File not on (or within 10 m of) the template grid: ", sources(r))
  }
  r <- shift(r, dx = shift[1], dy = shift[2])
  crs(r) <- crs(template)
  extend(crop(r, template), template)
}

build_monthly_ndvi_500m <- function(ndvi500_dir, template, out_dir, min_valid_frac = 0.25){

  if(is.character(template)) template <- rast(template)
  files <- list.files(ndvi500_dir, pattern = "^ndvi500_.*_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
  if(length(files) == 0) stop("No ndvi500_*.tif files in ", ndvi500_dir)

  months <- as.Date(paste0(sub("^.*_(\\d{4}-\\d{2})\\.tif$", "\\1", basename(files)), "-01"))
  if(anyDuplicated(months)) stop("More than one file for a month in ", ndvi500_dir)
  files <- files[order(months)]
  months <- sort(months)

  stacks <- lapply(c(mean = 1, sd = 2, valid = 3), function(b){
    # cells outside the template (sea, outside the study area) set to NA, as in build_monthly_ndvi()
    r <- mask(rast(lapply(files, function(f) on_template(rast(f, lyrs = b), template))), template[[1]])
    time(r) <- months
    names(r) <- format(months, "%Y-%m")
    r
  })

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stat_file <- function(s) file.path(out_dir, paste0(s, ".tif"))
  writeRaster(stacks$mean,  stat_file("raw_mean"), overwrite = TRUE)
  writeRaster(stacks$sd,    stat_file("raw_sd"), overwrite = TRUE)
  writeRaster(stacks$valid, stat_file("valid_fraction"), overwrite = TRUE)

  for(s in c("mean", "sd")){
    g <- fill_ndvi_gaps(stacks[[s]], stacks$valid, min_valid_frac)
    writeRaster(g$filled, stat_file(s), overwrite = TRUE)
    if(s == "mean") writeRaster(g$missing * 1, stat_file("filled"), overwrite = TRUE, datatype = "INT1U")
  }

  inside <- !is.na(values(template[[1]], mat = FALSE))
  filled <- values(rast(stat_file("filled")), mat = TRUE)[inside, , drop = FALSE]
  expected <- seq(min(months), max(months), by = "month")
  cat(sprintf("%s: %d months (%s to %s), %d without data, %.1f%% of land cell-months filled",
              out_dir, length(expected), format(min(months), "%Y-%m"), format(max(months), "%Y-%m"),
              sum(!expected %in% months), 100 * mean(filled)), fill = TRUE)
  invisible(out_dir)
}
