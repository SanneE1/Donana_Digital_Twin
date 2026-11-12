#
# Monthly flood map (0/1) on the large (lynx) template, for the rabbit K map (K = 0 where flood = 1)
#
#   flood_monthly_500m_large.tif   one layer per month (time = first of the month), 0/1 on land
#
# Sources, in order of preference:
#   1. inside the park (park_mask_500_large.tif): the LAST monthly flood maps the rabbit model was
#      calibrated with (LAST_monthly_floodmaps/<YYYY>_<MM>.tif, made by make_monthly_floodmaps()),
#      in the months that have one
#   2. elsewhere, and in park months without a LAST map: the Earth Engine water layers
#      (backcast/ndvi_pipeline/height_water_gee.py water): JRC Global Surface Water up to 2016
#      (water_fraction > 0), Sentinel-2 SCL water from 2017 (water_fraction > 0.05). Thresholds from
#      the comparison with LAST in the park (2004-2025): JRC > 0 kappa 0.63 (misses 40% of the
#      flooded cells: no data under clouds / Landsat 7 stripes), Sentinel-2 > 0.05 kappa 0.52.
#      Cells without a water observation count as not flooded.
#   3. every month: cells where rabbits cannot live (rabbit_no_habitat_500_large.tif: towns,
#      salines, water) are set to 1, so their K is 0.
#
# Run from the Donana_Digital_Twin project root after create_large_landcover_layers.R.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

data_dir <- file.path("data", "environmental_data")
template <- rast(file.path(data_dir, "template_raster_500_large.tif"))[[1]]
park     <- rast(file.path(data_dir, "park_mask_500_large.tif"))
no_rab   <- rast(file.path(data_dir, "rabbit_no_habitat_500_large.tif"))

# Settings that a calling script can set before sourcing this one (e.g. the 1985 backcast,
# backcast/Workflow_backcast_components.R):
#   flood_months     months to make
#   flood_last_dirs  folders with monthly LAST maps (<YYYY>_<MM>.tif); in a park month the first
#                    folder that has a map is used
#   flood_out_file   output file
if (!exists("flood_months"))    flood_months <- seq(as.Date("2004-01-01"), as.Date("2025-12-01"), by = "month")
if (!exists("flood_last_dirs")) flood_last_dirs <- file.path(data_dir, "LAST_monthly_floodmaps")
if (!exists("flood_out_file"))  flood_out_file <- file.path(data_dir, "flood_monthly_500m_large.tif")

months <- flood_months
jrc_until <- as.Date("2016-12-01")
thresholds <- c(jrc = 0, "sentinel-2" = 0.05)
gee_file <- function(src, m) file.path(data_dir, "water_large", src, sprintf("water500_%s_%s.tif", src, format(m, "%Y-%m")))

in_park <- values(park, mat = FALSE) %in% 1
land    <- !is.na(values(template, mat = FALSE))
blocked <- values(no_rab, mat = FALSE) %in% 1
source_count <- matrix(0L, length(months), 3, dimnames = list(format(months, "%Y-%m"), c("last", "gee", "none")))

layers <- lapply(seq_along(months), function(i){
  m <- months[i]
  flood <- rep(0, ncell(template))

  src <- if (m <= jrc_until) "jrc" else "sentinel-2"
  f <- gee_file(src, m)
  if (file.exists(f)) {
    g <- rast(f, lyrs = "water_fraction")
    g <- resample(g, template, method = "near")
    wf <- values(g, mat = FALSE)
    wf[wf < 0] <- NA                                   # NODATA -9999
    flood[!is.na(wf) & wf > thresholds[[src]]] <- 1
    source_count[i, "gee"] <<- 1L
  }

  lf <- file.path(flood_last_dirs, sprintf("%s_%s.tif", format(m, "%Y"), format(m, "%m")))
  lf <- lf[file.exists(lf)][1]
  if (!is.na(lf)) {
    l <- !is.na(values(resample(rast(lf), template, method = "near"), mat = FALSE))
    flood[in_park] <- as.numeric(l[in_park])
    source_count[i, "last"] <<- 1L
  }
  if (sum(source_count[i, 1:2]) == 0) source_count[i, "none"] <<- 1L

  flood[blocked] <- 1
  flood[!land] <- NA
  r <- template; values(r) <- flood
  r
})

flood <- rast(layers)
time(flood) <- months
names(flood) <- format(months, "%Y_%m")
writeRaster(flood, flood_out_file, overwrite = TRUE, datatype = "INT1U")

cat(sprintf("%d months: LAST in the park %d, Earth Engine %d, no water data %d\n", length(months),
            sum(source_count[, "last"]), sum(source_count[, "gee"]), sum(source_count[, "none"])))
if (any(source_count[, "none"] == 1)) cat("months without any water data (only the no-rabbit cells flooded):",
                                          rownames(source_count)[source_count[, "none"] == 1], "\n")
share <- global(flood, "mean", na.rm = TRUE)[, 1]
cat("flooded share of land cells per month: min", round(min(share), 3), "median", round(median(share), 3),
    "max", round(max(share), 3), "(no-rabbit cells alone", round(mean(blocked[land]), 3), ")\n")
