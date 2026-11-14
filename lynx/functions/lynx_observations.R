# ----------------------------------------------------------------------------
# Lynx observations and the statistics that compare the IBM with them
#
# Kept generic, so other presence or count sources can be used:
#   presence  a folder with one file per year, the year in the file name (e.g. 2005.shp,
#             presence_2005.gpkg, 2005.tif): polygons of the area where lynx were present, or a
#             raster (> 0 = present). Rasterised on the template (cells whose centre is inside).
#   counts    a CSV with columns year and breeding_females (other columns ignored)
# ----------------------------------------------------------------------------

# Presence per year: matrix cells x years (0/1, NA outside the template), colnames = years
read_lynx_presence <- function(dir, template, years = NULL){
  files <- list.files(dir, "\\.(shp|gpkg|geojson|tif)$", full.names = TRUE)
  yr <- as.integer(regmatches(basename(files), regexpr("(19|20)\\d{2}", basename(files))))
  keep <- !is.na(yr) & (is.null(years) | yr %in% years)
  files <- files[keep]; yr <- yr[keep]
  if (length(files) == 0) stop("No presence files (with a year in the name) in ", dir)
  land <- !is.na(terra::values(template[[1]], mat = FALSE))
  area <- terra::as.polygons(terra::ext(template), crs = terra::crs(template))
  out <- vapply(seq_along(files), function(i){
    f <- files[i]
    if (grepl("\\.tif$", f)) {
      r <- terra::resample(terra::rast(f), template, method = "max")
      v <- as.numeric(terra::values(r, mat = FALSE) > 0)
    } else {
      p <- terra::project(terra::vect(f), terra::crs(template))
      p <- terra::crop(p, area)
      v <- if (nrow(p) == 0) rep(0, terra::ncell(template))
           else as.numeric(terra::values(terra::rasterize(p, template), mat = FALSE) %in% 1)
    }
    v[is.na(v)] <- 0
    v[!land] <- NA
    v
  }, numeric(terra::ncell(template)))
  out <- out[, order(yr), drop = FALSE]
  colnames(out) <- sort(yr)
  out
}

read_lynx_counts <- function(file){
  d <- utils::read.csv(file)
  if (!all(c("year", "breeding_females") %in% names(d))) stop(file, " needs columns year and breeding_females")
  d[!is.na(d$breeding_females), c("year", "breeding_females")]
}

# Aggregates cell values (matrix cells x years, template order) to coarser cells of fact x fact
# template cells: 1 where any cell is 1. Returns a matrix coarse cells x years.
coarse_presence <- function(m, template, fact){
  r <- terra::rast(template[[1]], nlyrs = ncol(m))
  terra::values(r) <- m
  a <- terra::aggregate(r, fact, max, na.rm = TRUE)
  terra::values(a)
}

# Statistics of one simulation against the observations, over `years`:
#   d_count     RMSE of log(1 + territorial females) - log(1 + observed breeding females)
#               With late_from (a year), split in d_count_early (years before late_from) and
#               d_count_late (late_from on), so that a run is judged on the later level of the
#               population separately (the plateau of the Doñana population from about 2012):
#               otherwise a run that is too low early and too high late can fit as well as one
#               that follows the counts.
#   d_presence  1 - mean Dice overlap (2 |S n O| / (|S| + |O|)) of simulated territories and
#               observed presence, both on coarse cells (fact x fact template cells, any = 1);
#               1 when the simulation has no territories
# plus, for evaluation, the per-year values.
lynx_fit_stats <- function(sim_terr, sim_occ, obs_counts, obs_presence, template, years, fact = 4,
                           late_from = NULL){
  cy <- intersect(years, obs_counts$year)
  sim_f <- sim_terr$females_territorial[match(cy, sim_terr$year)]
  sim_f[is.na(sim_f)] <- 0
  obs_f <- obs_counts$breeding_females[match(cy, obs_counts$year)]
  rmse <- function(k) if (any(k)) sqrt(mean((log1p(sim_f[k]) - log1p(obs_f[k]))^2)) else NA
  d_count <- if (is.null(late_from)) c(d_count = rmse(rep(TRUE, length(cy)))) else
    c(d_count_early = rmse(cy < late_from), d_count_late = rmse(cy >= late_from))

  py <- intersect(intersect(years, as.integer(colnames(obs_presence))), as.integer(colnames(sim_occ)))
  dice <- numeric(0)
  if (length(py) > 0) {
    s <- coarse_presence(sim_occ[, as.character(py), drop = FALSE], template, fact)
    o <- coarse_presence(obs_presence[, as.character(py), drop = FALSE], template, fact)
    dice <- vapply(seq_along(py), function(j){
      sj <- s[, j] %in% 1; oj <- o[, j] %in% 1
      if (sum(sj) + sum(oj) == 0) return(1)
      2 * sum(sj & oj) / (sum(sj) + sum(oj))
    }, 0)
    names(dice) <- py
  }
  list(stats = c(d_count, d_presence = if (length(dice)) 1 - mean(dice) else NA),
       counts = data.frame(year = cy, sim = sim_f, obs = obs_f),
       dice = dice)
}
