#
# OpenStreetMap roads of the large (lynx) template area
#
# Reads the Spain OpenStreetMap extract (Geofabrik .osm.pbf, in the PVA project) once and keeps the
# public roads inside the large template (+ 2 km), so the lynx maps do not need the 1.4 GB file:
#   data/environmental_data/roads/osm_roads_large.gpkg   lines with the OSM `highway` class
#
# Run from the Donana_Digital_Twin project root after create_large_template.R.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)
library(sf)

data_dir <- file.path("data", "environmental_data")
template <- rast(file.path(data_dir, "template_raster_500_large.tif"))
pbf_file <- "C:/Users/z1512834z/Documents/01_Natura_Connect/04_Lynx_PVA/data/original_data/OpenStreetMapRoads/spain-260712.osm.pbf"
out_file <- file.path(data_dir, "roads", "osm_roads_large.gpkg")

road_classes <- c("motorway", "motorway_link", "trunk", "trunk_link",
                  "primary", "primary_link", "secondary", "secondary_link",
                  "tertiary", "tertiary_link", "unclassified")

area <- st_as_sfc(st_bbox(c(xmin = xmin(template) - 2000, ymin = ymin(template) - 2000,
                            xmax = xmax(template) + 2000, ymax = ymax(template) + 2000),
                          crs = st_crs(crs(template))))
area_ll <- st_transform(area, 4326)

query <- paste0("SELECT highway FROM lines WHERE highway IN (",
                paste0("'", road_classes, "'", collapse = ", "), ")")
roads <- st_read(pbf_file, layer = "lines", query = query, wkt_filter = st_as_text(area_ll), quiet = TRUE)
roads <- st_intersection(st_transform(roads, st_crs(area)), area)
roads <- roads[st_geometry_type(roads) %in% c("LINESTRING", "MULTILINESTRING"), ]

dir.create(dirname(out_file), showWarnings = FALSE)
st_write(roads, out_file, delete_dsn = TRUE, quiet = TRUE)
cat(nrow(roads), "road lines written to", out_file, "\n")
print(round(tapply(as.numeric(st_length(roads)), roads$highway, sum) / 1000))
