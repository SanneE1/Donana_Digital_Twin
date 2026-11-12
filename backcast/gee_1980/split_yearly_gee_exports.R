#
# Split the yearly multi-band exports of gee_1980/01, 02 and 05 (*.js) into the per-month files
# that the existing R code reads:
#   ndvi500_<sensor>_<YYYY>.tif  ->  ndvi500_<sensor>_<YYYY-MM>.tif  (bands mean, sd, valid_fraction;
#                                    read by build_monthly_ndvi_500m())
#   water500_jrc_<YYYY>.tif      ->  water500_jrc_<YYYY-MM>.tif      (bands water_fraction,
#                                    valid_fraction; read by create_large_flood_layer.R)
# Months without any observation (valid_fraction 0 everywhere) are not written, as the Python
# exporters do not export months without scenes (they become gaps, filled later).
#
#   Rscript backcast/gee_1980/split_yearly_gee_exports.R <folder with yearly files> <output folder>

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "^(ndvi500|water500)_.+_\\d{4}(-\\d+-\\d+)?\\.tif$", full.names = TRUE)
for (f in files) {
  r <- rast(f)
  # GEE splits large exports into tiles (name-0000000000-0000000000.tif); the template is small
  # enough for one file, but merge if needed before running this.
  stem <- sub("_(\\d{4})(-\\d+-\\d+)?\\.tif$", "", basename(f))
  ym <- unique(sub("^(\\d{4})_(\\d{2})_.*$", "\\1-\\2", names(r)))
  for (m in ym) {
    lyr <- r[[startsWith(names(r), paste0(sub("-", "_", m), "_"))]]
    names(lyr) <- sub("^\\d{4}_\\d{2}_", "", names(lyr))
    if (all(values(lyr[["valid_fraction"]], mat = FALSE) %in% c(0, NA))) next
    lyr <- lyr[[if (startsWith(stem, "ndvi500")) c("mean", "sd", "valid_fraction")
                else c("water_fraction", "valid_fraction")]]
    writeRaster(lyr, file.path(out_dir, sprintf("%s_%s.tif", stem, m)), overwrite = TRUE,
                NAflag = -9999)
  }
  cat(basename(f), "->", length(ym), "months checked\n")
}
