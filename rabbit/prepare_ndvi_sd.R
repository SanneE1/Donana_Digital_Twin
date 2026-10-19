
# One-off preparation: monthly 500 m NDVI mean and standard-deviation rasters
#
# The rabbit model (create_K_map_gaussian) needs a monthly NDVI mean and SD covariate. They are
# derived exactly as in Workflow_calibration.R (10 m -> 500 m aggregation with fun = "mean" / "sd",
# then monthly max), but run one file at a time because the raw NDVI folder is ~19 GB.
# NB: data/environmental_data/monthly_ndvi.tif is NOT the same as this mean (it is a different
# monthly composite), so do not use it with the calibrated rabbit parameters.
#
# Output (small, can be shipped with the shiny app):
#   data/environmental_data/monthly_ndvi_mean.tif
#   data/environmental_data/monthly_ndvi_sd.tif
#
# Run from the Donana_Digital_Twin project root. Existing outputs are kept; delete them to recompute.

library(terra)
library(ncdf4)
library(lubridate)

source(file.path("rabbit", "functions", "miscellaneous.R"))

template   <- rast(file.path("data", "environmental_data", "template_raster_500.tif"))
ndvi_files <- list.files(file.path("data", "environmental_data", "NDVI"),
                         pattern = "\\.nc$", full.names = TRUE)

aggregate_ndvi <- function(fun, out_file) {
  if (file.exists(out_file)) {
    cat(out_file, "already exists, skipping\n")
    return(invisible(NULL))
  }

  r_list <- list()
  for (f in ndvi_files) {
    cat("Processing", basename(f), "(", fun, ")\n")
    r_list[[f]] <- load_and_format_rasters(f, template,
                                           aggregate = TRUE, fun = fun,
                                           temporal_agg = "month", temporal_fun = "max")
  }

  # A month can be split over two files; keep the maximum, as the monthly aggregation does
  r <- rast(r_list)
  tt <- as.Date(time(r))
  if (anyDuplicated(tt)) {
    r <- terra::tapp(r, index = tt, fun = "max", na.rm = TRUE)
    time(r) <- sort(unique(tt))
  }
  r <- r[[order(time(r))]]

  writeRaster(r, out_file, overwrite = TRUE)
  cat("Written", out_file, "with", nlyr(r), "layers:",
      format(min(time(r))), "to", format(max(time(r))), "\n")
}

aggregate_ndvi("sd",   file.path("data", "environmental_data", "monthly_ndvi_sd.tif"))
aggregate_ndvi("mean", file.path("data", "environmental_data", "monthly_ndvi_mean.tif"))
