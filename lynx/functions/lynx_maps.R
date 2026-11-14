# ----------------------------------------------------------------------------
# Yearly map inputs of the lynx IBM (lynx/Pascal_program) on the large template
#
# The IBM reads in June of every year <folder>/Lynx_HabitatMap_<year>.txt (0 = barrier,
# 1 = matrix/open, 2 = dispersal habitat), Lynx_BreedingMap_<year>.txt and Lynx_PreyMap_<year>.txt
# (0/1). Format: a first line "ncol nrow", then one line per map row from north to south; values
# < 0 (no data) are read as 0.
#
# The maps of year y describe the twelve months before that June: June of y - 1 to May of y.
# ----------------------------------------------------------------------------

# Months (Date, first of the month) whose data make the maps of `year`
lynx_year_months <- function(year){
  seq(as.Date(sprintf("%04d-06-01", year - 1)), by = "month", length.out = 12)
}

# Layers of a monthly stack (time = first of the month) for the months of `year`; NULL if none
select_year <- function(r, year){
  i <- which(format(as.Date(terra::time(r)), "%Y-%m") %in% format(lynx_year_months(year), "%Y-%m"))
  if (length(i) == 0) NULL else r[[i]]
}

# Habitat map from NDVI and vegetation height, per year. Barrier (towns) and road cells are the same
# as in the CORINE map (environmental_data/create_large_landcover_layers.R), so the two maps differ
# only in what counts as cover.
#
#   ndvi_mean, ndvi_sd  monthly 500 m NDVI mean and within-cell SD (the rabbit model's covariates,
#                       gap-filled, e.g. NDVI_large/harmonised/mean.tif and sd.tif)
#   barrier, roads      0/1 layers: towns (barrier, 0) and cells crossed by roads (matrix, 1, never
#                       habitat): lynx_barrier_500_large.tif, lynx_roads_500_large.tif
#   height              static mean vegetation height (m), or NULL
#   p                   list of thresholds (calibrated):
#     t_cover  habitat where the year's maximum NDVI >= t_cover ...
#     sd_lo, sd_hi  ... and the year's mean within-cell NDVI SD is in [sd_lo, sd_hi],
#     t_h      or where height >= t_h (whatever the NDVI); Inf to switch off
#   Everything else (open land, crops, marsh, water) is matrix.
#
# Returns a SpatRaster with one layer per year (names = years), NA outside the template.
lynx_habitat_ndvi <- function(ndvi_mean, ndvi_sd, years, p, template, barrier, roads, height = NULL){
  out <- lapply(years, function(y){
    m <- select_year(ndvi_mean, y); s <- select_year(ndvi_sd, y)
    if (is.null(m) || is.null(s)) stop("No NDVI for the lynx year ", y, " (June ", y - 1, " to May ", y, ")")
    nmax <- max(m, na.rm = TRUE)
    smean <- mean(s, na.rm = TRUE)

    cover <- nmax >= p$t_cover & smean >= p$sd_lo & smean <= p$sd_hi
    if (!is.null(height) && is.finite(p$t_h)) cover <- cover | (height >= p$t_h)
    hab <- terra::ifel(barrier == 1, 0, terra::ifel(cover & roads == 0, 2, 1))
    terra::mask(hab, template[[1]])
  })
  out <- terra::rast(out)
  names(out) <- years
  out
}

# Prey maps from monthly rabbit density: 1 where the density is >= threshold in at least n_months
# of the twelve months of the lynx year (June y - 1 to May y), as in the PVA project
# (Create_breeding_maps.R). `density` is a monthly stack (time = first of the month), or a
# matrix cells x months with `months` (Date) given, which is much faster inside a calibration.
# Returns a matrix cells x years (0/1, NA outside the template) with colnames = years.
lynx_prey_matrix <- function(density, years, threshold, n_months, months = NULL){
  if (inherits(density, "SpatRaster")) {
    months <- as.Date(terra::time(density))
    density <- terra::values(density)
  }
  ym <- format(months, "%Y-%m")
  out <- vapply(years, function(y){
    i <- which(ym %in% format(lynx_year_months(y), "%Y-%m"))
    if (length(i) < 12) stop("Rabbit density lacks months of the lynx year ", y, " (June ", y - 1, " to May ", y, ")")
    as.numeric(rowSums(density[, i, drop = FALSE] >= threshold) >= n_months)
  }, numeric(nrow(density)))
  out[is.na(density[, 1]), ] <- NA
  colnames(out) <- years
  out
}

# Map from a vector of cell values on the template
cells_to_map <- function(values, template){
  r <- terra::rast(template[[1]])
  terra::values(r) <- values
  r
}

# Breeding map: every land cell (territories then depend on habitat = 2 and the prey map only)
lynx_breeding_all <- function(template){
  terra::mask(terra::init(terra::rast(template[[1]]), 1), template[[1]])
}

# Population ID map: survival in the IBM depends on the natal population ID. 4 (Doñana, survival
# column 2) inside the park, 0 (survival column 0) elsewhere, so lynx born inside and outside the
# park can have different survival, as in Revilla & Wiegand 2008.
lynx_pops_map <- function(template, park){
  r <- terra::ifel(!is.na(park) & park == 1, 4, 0)
  terra::mask(r, template[[1]])
}

# Writes one map in the IBM text format
write_lynx_map <- function(r, file){
  v <- terra::values(r[[1]], mat = FALSE)
  v[is.na(v)] <- -9999
  m <- matrix(as.integer(round(v)), nrow = terra::nrow(r), byrow = TRUE)
  con <- file(file, "w"); on.exit(close(con))
  writeLines(paste(terra::ncol(r), terra::nrow(r)), con)
  utils::write.table(m, con, row.names = FALSE, col.names = FALSE)
  invisible(file)
}

# Writes the layers of `maps` (names = years) as <dir>/<prefix>_<year>.txt, e.g. prefix
# "Lynx_HabitatMap". A static map (one layer) is written for every year in `years`.
write_lynx_map_series <- function(maps, dir, prefix, years = as.integer(names(maps))){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (y in years) {
    r <- if (terra::nlyr(maps) == 1) maps else maps[[as.character(y)]]
    write_lynx_map(r, file.path(dir, sprintf("%s_%d.txt", prefix, y)))
  }
  invisible(dir)
}
