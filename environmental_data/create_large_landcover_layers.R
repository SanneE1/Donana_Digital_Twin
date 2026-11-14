#
# Static land-cover layers on the large (lynx) template, from CORINE Land Cover 2018 and
# OpenStreetMap roads
#
#   corine2018_500_large.tif        CORINE class index (1-44, 48) per cell, majority of the 100 m
#                                   CORINE cells
#   lynx_barrier_500_large.tif      1 = town (barrier for lynx in every habitat map): CORINE built-up
#                                   classes (1, 2 urban fabric, 3 industrial/commercial, 5 ports,
#                                   6 airports) are the majority of the cell; 0 elsewhere
#   lynx_roads_500_large.tif        1 = crossed by a road of `road_classes` (OSM); such cells are
#                                   open matrix, never habitat, so roads are crossed but never inside
#                                   a territory; 0 elsewhere
#   lynx_habitat_corine_500_large.tif
#                                   lynx habitat map for the IBM: 0 = barrier, 1 = matrix/open,
#                                   2 = dispersal habitat. Habitat = the Revilla et al. 2015 cover
#                                   classes (forests 23-25, moors, sclerophyllous vegetation,
#                                   woodland-shrub 27-29; option 1 of the PVA project's
#                                   reclassify_corine.R), majority of the 100 m cells. Unlike Revilla
#                                   2015, marshes, salines and water are matrix (lynx move along and
#                                   across the marisma and can cross water) and only towns are
#                                   barrier; road cells are matrix.
#   rabbit_no_habitat_500_large.tif 1 where rabbits cannot live (K = 0): artificial surfaces except
#                                   green urban / sport areas (CORINE 1-9), salines and intertidal
#                                   flats (38, 39), water (40-44) and no data (48); 0 elsewhere.
#                                   Outside the park the rabbit K map is extrapolated from NDVI and
#                                   vegetation height, which would otherwise give rabbits in towns
#                                   and reservoirs.
#
# CORINE layers use the 500 m majority, so features smaller than half a cell (villages, narrow
# rivers) are not resolved. Run from the Donana_Digital_Twin project root after
# create_large_template.R and extract_osm_roads_large.R.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
library(terra)

data_dir <- file.path("data", "environmental_data")
template <- rast(file.path(data_dir, "template_raster_500_large.tif"))

# CORINE of the PVA project (not in this repository)
corine_file <- "C:/Users/z1512834z/Documents/01_Natura_Connect/04_Lynx_PVA/data/original_data/U2018_CLC2018_V2020_20u1.tif"
roads_file  <- file.path(data_dir, "roads", "osm_roads_large.gpkg")

town_classes  <- c(1, 2, 3, 5, 6)
cover_classes <- c(23:25, 27:29)
# OSM roads that make a cell open matrix. Share of CORINE habitat cells in the 2002-2017 lynx
# presence area that this removes: motorway/trunk 0.7%, + primary 1.0%, + secondary 4.0%,
# + tertiary 7.3%, + unclassified 13.6%.
road_classes <- c("motorway", "motorway_link", "trunk", "trunk_link",
                  "primary", "primary_link", "secondary", "secondary_link")

# CORINE (EPSG:3035, 100 m) to the template CRS at 100 m (nearest), then the majority per 500 m cell
corine <- rast(corine_file)
area   <- project(as.polygons(ext(template), crs = crs(template)), crs(corine))
grid100 <- rast(ext(template), resolution = 100, crs = crs(template))
c100   <- project(crop(corine, buffer(area, 1000)), grid100, method = "near")
c500   <- mask(resample(c100, template, method = "mode"), template[[1]])
names(c500) <- "corine2018"
writeRaster(c500, file.path(data_dir, "corine2018_500_large.tif"), overwrite = TRUE, datatype = "INT1U")

codes <- c(1:44, 48)
town <- mask(classify(c500, cbind(codes, as.integer(codes %in% town_classes))), template[[1]])
names(town) <- "lynx_barrier"
writeRaster(town, file.path(data_dir, "lynx_barrier_500_large.tif"), overwrite = TRUE, datatype = "INT1U")

roads <- vect(roads_file)
roads <- roads[roads$highway %in% road_classes]
road_cells <- mask(rasterize(roads, template, field = 1, background = 0, touches = TRUE), template[[1]])
names(road_cells) <- "lynx_roads"
writeRaster(road_cells, file.path(data_dir, "lynx_roads_500_large.tif"), overwrite = TRUE, datatype = "INT1U")

# Cover: majority of the 100 m cells cover vs not cover
cover <- resample(classify(c100, cbind(codes, as.integer(codes %in% cover_classes))), template, method = "mode")
habitat <- ifel(town == 1, 0, ifel(cover == 1 & road_cells == 0, 2, 1))
habitat <- mask(habitat, template[[1]])
names(habitat) <- "lynx_habitat"
writeRaster(habitat, file.path(data_dir, "lynx_habitat_corine_500_large.tif"), overwrite = TRUE, datatype = "INT1U")

no_rabbit <- mask(classify(c500, cbind(codes, as.integer(codes %in% c(1:9, 38:44, 48)))), template[[1]])
names(no_rabbit) <- "rabbit_no_habitat"
writeRaster(no_rabbit, file.path(data_dir, "rabbit_no_habitat_500_large.tif"), overwrite = TRUE, datatype = "INT1U")

land <- !is.na(values(template[[1]], mat = FALSE))
h <- values(habitat, mat = FALSE)[land]
cat(sprintf("lynx habitat (land cells): barrier %.3f, matrix %.3f, habitat %.3f; road cells %d\n",
            mean(h == 0), mean(h == 1), mean(h == 2), sum(values(road_cells, mat = FALSE)[land] == 1)))
cat(sprintf("rabbit no-habitat cells: %d of %d land cells (%.1f%%)\n",
            sum(values(no_rabbit, mat = FALSE)[land] == 1), sum(land),
            100 * mean(values(no_rabbit, mat = FALSE)[land] == 1)))
