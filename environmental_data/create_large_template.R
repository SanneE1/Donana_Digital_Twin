#
# Template raster for the larger Doñana area (lynx population), 500 m
#
# The rabbit model's template (template_raster_500.tif) covers the National Park only. The Doñana
# lynx population ranges well beyond it, so the lynx model (and the NDVI and other drivers it needs)
# uses a larger grid:
#   - extent: all Doñana lynx presence cells of the annual distribution maps (1 km, 2002-2017) and
#     the buffered 2022 Doñana population, plus `buffer_m`
#   - same CRS, resolution and cell alignment as template_raster_500.tif, so every park cell is
#     exactly one cell of the large grid
#   - value 1 on land, NA at sea (cells that are mostly CORINE class 44, "sea and ocean", or the
#     CORINE NODATA of the open sea)
#
# Writes to data/environmental_data:
#   template_raster_500_large.tif     the template
#   park_mask_500_large.tif           1 for the cells of template_raster_500.tif, NA elsewhere
#   Donana_large_template_extent/     extent of the template as a shapefile (download area)
#
# Run from the Donana_Digital_Twin project root.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

data_dir      <- file.path("data", "environmental_data")
park_template <- rast(file.path(data_dir, "template_raster_500.tif"))

# Lynx data of the PVA project (not in this repository)
pva_data       <- "C:/Users/z1512834z/Documents/01_Natura_Connect/04_Lynx_PVA/data"
distribution   <- file.path(pva_data, "original_data", "Annual_distribution_Alejandro")
distr_years    <- 2002:2017
pops_2022_file <- file.path(pva_data, "GIS_maps", "Lynx_populations_2022_buffered_vector.shp")
corine_file    <- file.path(pva_data, "original_data", "U2018_CLC2018_V2020_20u1.tif")

search_m <- 60000   # only lynx data within this distance of the park grid count as Doñana
buffer_m <- 5000    # margin around the Doñana lynx data
sea_class <- 44     # CORINE "sea and ocean"

out_template <- file.path(data_dir, "template_raster_500_large.tif")
out_park     <- file.path(data_dir, "park_mask_500_large.tif")
out_extent   <- file.path(data_dir, "Donana_large_template_extent", "large_template_extent.shp")

#-------------------------------------------------------------------------------
# Extent of the Doñana lynx population
#-------------------------------------------------------------------------------

near <- buffer(as.polygons(ext(park_template), crs = crs(park_template)), search_m)

doñana_part <- function(v){
  v <- project(v, crs(park_template))
  v[relate(v, near, "intersects")[, 1]]
}

presence <- lapply(distr_years, function(y){
  doñana_part(vect(file.path(distribution, paste0(y, "_Iberian lynx distribution.shp"))))
})
presence <- presence[sapply(presence, nrow) > 0]
pops_2022 <- doñana_part(vect(pops_2022_file))

lynx_ext <- Reduce(union, c(lapply(presence, ext), list(ext(pops_2022))))
cat("Doñana lynx data extent:", as.vector(lynx_ext), fill = TRUE)

#-------------------------------------------------------------------------------
# Grid aligned with the park template
#-------------------------------------------------------------------------------

e  <- ext(park_template)
rs <- res(park_template)[1]
grow <- function(edge, target, direction) edge + direction * ceiling(direction * (target - edge) / rs) * rs

large_ext <- ext(grow(e$xmin, lynx_ext$xmin - buffer_m, -1),
                 grow(e$xmax, lynx_ext$xmax + buffer_m,  1),
                 grow(e$ymin, lynx_ext$ymin - buffer_m, -1),
                 grow(e$ymax, lynx_ext$ymax + buffer_m,  1))
# the park grid is never cut off
large_ext <- union(large_ext, e)

template <- rast(large_ext, resolution = rs, crs = crs(park_template))

#-------------------------------------------------------------------------------
# Land / sea
#-------------------------------------------------------------------------------

corine <- crop(rast(corine_file), project(as.polygons(ext(template), crs = crs(template)), crs(rast(corine_file))), snap = "out")
# CORINE only maps a strip of coastal sea as class 44; the open sea further out is NODATA (48 or NA)
sea <- ifel(is.na(corine) | corine == sea_class | corine == 48, 1, 0)
sea_frac <- project(sea, template, method = "average")
values(template) <- ifelse(values(sea_frac) > 0.5, NA, 1)

#-------------------------------------------------------------------------------
# Checks and output
#-------------------------------------------------------------------------------

park_mask <- extend(park_template, template)
park_mask[!is.na(park_mask)] <- 1
stopifnot(compareGeom(park_mask, template, stopOnError = FALSE),
          sum(!is.na(values(park_mask))) == sum(!is.na(values(park_template))))

writeRaster(template, out_template, overwrite = TRUE, NAflag = -9999)
writeRaster(park_mask, out_park, overwrite = TRUE, NAflag = -9999)
dir.create(dirname(out_extent), showWarnings = FALSE)
writeVector(as.polygons(ext(template), crs = crs(template)), out_extent, overwrite = TRUE)

cat(sprintf("Large template: %d x %d cells (%.0f x %.0f km), %d land cells (%.0f km2); park grid: %d x %d cells",
            ncol(template), nrow(template), ncol(template) * rs / 1000, nrow(template) * rs / 1000,
            sum(!is.na(values(template))), sum(!is.na(values(template))) * rs^2 / 1e6,
            ncol(park_template), nrow(park_template)), fill = TRUE)
cat("Extent:", as.vector(ext(template)), fill = TRUE)
