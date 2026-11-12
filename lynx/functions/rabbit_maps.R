# ----------------------------------------------------------------------------
# Rabbit density maps for the lynx model
#
# The rabbit input can be
#   - one file with one layer per month (time = first of the month), e.g. the mean of a rabbit
#     hindcast, or a directory holding only that file;
#   - a directory in the layout of the rabbit calibration's pop_maps/ (its in-sample output
#     historic_in_sample_simulation/, or pop_maps/ itself) or of the spatial projection's
#     map_series/ (census-scaled maps):
#         runs/run_<k>/<map>_<yyyy>_<m>.tif      one map per month and run
#         summary/<map>_<yyyy>_<m>__mean.tif     the mean over the runs
#
# mode "mean": one density for all simulations: the file, the summary means, or the mean over
#   the runs when there is no summary.
# mode "runs": the maps of every run (not with a single file). Each lynx simulation uses one
#   run, drawn at random with its seed (lynx_simulate()), so the uncertainty of the rabbit maps
#   enters the lynx calibration.
# ----------------------------------------------------------------------------

# Months (first of the month, Date) from map names ending in _<yyyy>_<m><suffix>; NA otherwise
map_name_months <- function(files, suffix = ""){
  base <- tools::file_path_sans_ext(basename(files))
  ym <- regmatches(base, regexec(paste0("_(\\d{4})_(\\d{1,2})", suffix, "$"), base))
  as.Date(vapply(ym, function(x) if (length(x) == 3) sprintf("%s-%02d-01", x[2], as.integer(x[3])) else NA_character_,
                 character(1)))
}

# Where the maps are: list(file) or list(runs, summary) (directories, may not exist)
rabbit_map_layout <- function(path){
  if (!dir.exists(path)) {
    if (!file.exists(path)) stop("Rabbit maps not found: ", path, call. = FALSE)
    return(list(file = path))
  }
  for (d in c(file.path(path, "pop_maps"), path)) {
    if (dir.exists(file.path(d, "runs")) || dir.exists(file.path(d, "summary"))) {
      return(list(runs = file.path(d, "runs"), summary = file.path(d, "summary")))
    }
  }
  tifs <- list.files(path, pattern = "[.]tif$", full.names = TRUE)
  if (length(tifs) == 1) return(list(file = tifs))
  stop("No rabbit maps in ", path, ": expected one monthly .tif file, or runs/ and/or summary/ ",
       "(the layout of the rabbit calibration's pop_maps/ or the spatial projection's map_series/)",
       call. = FALSE)
}

# Values of the maps of `months` (one file per month) on the land cells: matrix land cells x months
read_month_maps <- function(files, file_months, months, template, land, what){
  miss <- months[!format(months, "%Y-%m") %in% format(file_months, "%Y-%m")]
  if (length(miss)) {
    stop(what, " has no map for ", length(miss), " month(s) the lynx years need (",
         paste(format(utils::head(miss, 5), "%Y-%m"), collapse = ", "), if (length(miss) > 5) ", ...",
         "). The rabbit calibration writes every month since its maps are saved per run; ",
         "older outputs only have March, June and September", call. = FALSE)
  }
  f <- files[match(format(months, "%Y-%m"), format(file_months, "%Y-%m"))]
  r <- terra::rast(f)
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) stop(what, ": maps not on the template grid", call. = FALSE)
  v <- terra::values(r)[land, , drop = FALSE]
  if (all(v == round(v), na.rm = TRUE)) storage.mode(v) <- "integer"   # half the memory
  dimnames(v) <- NULL
  v
}

# path: file or directory (see above); template: path or SpatRaster; months: the months needed
# (Date, first of the month), NULL = all months of the maps.
# Returns an object of class rabbit_maps: list(mode, months, land (logical per cell), and
#   density (cells x months, NA off land) for mode mean, or runs (list of land cells x months
#   matrices) and run_names for mode runs, and source (a description)).
read_rabbit_maps <- function(path, template, mode = c("mean", "runs"), months = NULL){
  mode <- match.arg(mode)
  if (is.character(template)) template <- terra::rast(template)
  template <- template[[1]]
  land <- !is.na(terra::values(template, mat = FALSE))
  lay <- rabbit_map_layout(path)
  out <- structure(list(mode = mode, land = land), class = "rabbit_maps")
  full <- function(v){
    m <- matrix(NA_real_, length(land), ncol(v))
    m[land, ] <- v
    m
  }

  if (!is.null(lay$file)) {
    if (mode == "runs") stop("rabbitMaps = runs needs the maps per run (runs/), not one file", call. = FALSE)
    r <- terra::rast(lay$file)
    file_months <- as.Date(terra::time(r))
    if (any(is.na(file_months))) stop("Rabbit density file without a time per layer: ", lay$file, call. = FALSE)
    if (is.null(months)) months <- file_months
    miss <- months[!format(months, "%Y-%m") %in% format(file_months, "%Y-%m")]
    if (length(miss)) stop("The rabbit density lacks ", length(miss), " month(s) the lynx years need (from ",
                           format(min(miss), "%Y-%m"), ")", call. = FALSE)
    r <- r[[match(format(months, "%Y-%m"), format(file_months, "%Y-%m"))]]
    if (!terra::compareGeom(r, template, stopOnError = FALSE)) stop("Rabbit density not on the template grid", call. = FALSE)
    out$density <- terra::values(r)
    out$density[!land, ] <- NA
    out$months <- months
    out$source <- basename(lay$file)
    return(out)
  }

  run_dirs <- if (dir.exists(lay$runs)) list.dirs(lay$runs, recursive = FALSE) else character(0)
  run_files <- lapply(run_dirs, list.files, pattern = "[.]tif$", full.names = TRUE)
  keep <- lengths(run_files) > 0
  run_dirs <- run_dirs[keep]; run_files <- run_files[keep]
  run_months <- lapply(run_files, map_name_months)
  if (is.null(months) && length(run_files)) months <- sort(Reduce(intersect, run_months))
  if (!is.null(months)) months <- as.Date(months)

  if (mode == "mean") {
    sfiles <- if (dir.exists(lay$summary)) list.files(lay$summary, pattern = "__mean[.]tif$", full.names = TRUE) else character(0)
    if (length(sfiles)) {
      smonths <- map_name_months(sfiles, "__mean")
      if (is.null(months)) months <- sort(smonths[!is.na(smonths)])
      out$density <- full(read_month_maps(sfiles, smonths, months, template, land, "The rabbit summary (mean)"))
      out$source <- "summary means"
    } else {
      if (!length(run_files)) stop("No rabbit maps in ", lay$runs, " or ", lay$summary, call. = FALSE)
      s <- 0
      for (i in seq_along(run_files)) {
        s <- s + read_month_maps(run_files[[i]], run_months[[i]], months, template, land, basename(run_dirs[i]))
      }
      out$density <- full(s / length(run_files))
      out$source <- sprintf("mean of %d runs", length(run_files))
    }
  } else {
    if (!length(run_files)) stop("rabbitMaps = runs: no maps per run in ", lay$runs, call. = FALSE)
    out$runs <- lapply(seq_along(run_files), function(i){
      read_month_maps(run_files[[i]], run_months[[i]], months, template, land, basename(run_dirs[i]))
    })
    out$run_names <- basename(run_dirs)
    out$source <- sprintf("%d runs", length(run_files))
  }
  out$months <- months
  out
}

# The months of the lynx years start_year to end_year (June of the year before to May)
lynx_years_months <- function(start_year, end_year){
  seq(as.Date(sprintf("%04d-06-01", start_year - 1)), as.Date(sprintf("%04d-05-01", end_year)), by = "month")
}

# Prior range of prey_threshold on the scale of the rabbit maps of `inp` (prepare_lynx_inputs()).
# Rabbit maps differ in level (simulated in-sample maps, census-scaled maps of spatial-projection,
# a hindcast), so a fixed prior would cover different parts of the density distribution. The
# range is c(1/40, 1) x the 99.5% quantile of the monthly density of the land cells in the lynx
# years (for runs: of all runs together). On the rabbit hindcast this gives the earlier fixed
# prior, 0.5-20 rabbits per cell (20 is exceeded in about 0.4% of the cell-months).
prey_threshold_range <- function(inp, q = 0.995, ratio = 40){
  thin <- function(v, n) v[seq(1, length(v), by = max(1, floor(length(v) / n)))]
  v <- if (!is.null(inp$density_runs)) {
    unlist(lapply(inp$density_runs, function(m) thin(as.vector(m), 1e6 / length(inp$density_runs))))
  } else thin(as.vector(inp$density[inp$land, ]), 1e6)
  hi <- as.numeric(stats::quantile(v, q, na.rm = TRUE))
  if (!is.finite(hi) || hi <= 0) stop("The rabbit maps are (almost) all zero: no prior range for prey_threshold", call. = FALSE)
  signif(c(hi / ratio, hi), 3)
}
