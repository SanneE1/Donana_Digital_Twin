# ----------------------------------------------------------------------------
# One lynx IBM simulation for a parameter set, on inputs prepared once by prepare_lynx_inputs()
#
# Calibrated parameters (lynx_priors()):
#   Tsize           female territory size (500 m cells; males x male_T_multiplier), rounded
#   surv_out        annual survival from age 2 of lynx born outside the park
#   surv_disp_rho   extra daily mortality of dispersers per share of steps in open matrix
#                   (narrow prior around the rates file's 5.8)
#   prey_threshold  rabbit density (rabbits per 500 m cell, monthly mean of the rabbit
#                   hindcast/forecast) a territory cell needs ...
#   prey_months     ... in at least this many months of the lynx year (June-May), rounded
#   start_mult      start lynx = start_mult x the observed breeding females of start_year: the
#                   start lynx settle as adults, while the real population also had subadults
#                   and young; with 2x the simulated population stays near its start size
# Fixed at the rates file (not calibrated; the counts cannot separate them from surv_out and
# prey): rep_prob (0.57) and the age-specific survival of lynx born inside the park (Doñana
# column of surv_<age>). A posterior that does contain rep_prob or surv_in (calibrated before
# 2026-11-09) still uses its own values.
# and for the NDVI habitat map (map_type "ndvi", lynx_habitat_ndvi()):
#   t_cover         habitat where the year's maximum NDVI >= t_cover ...
#   sd_hi           ... and the year's mean within-cell NDVI SD <= sd_hi,
#   t_h             or where vegetation height >= t_h (m)
# ----------------------------------------------------------------------------

# prey_range: prior range of prey_threshold, on the scale of the rabbit maps used
# (prey_threshold_range(), rabbit_maps.R); 0.5-20 rabbits per cell fits the rabbit hindcast
lynx_priors <- function(map_type = c("corine", "ndvi"), prey_range = c(0.5, 20)){
  map_type <- match.arg(map_type)
  p <- list(
    Tsize          = c("unif", 4, 24.49),   # 1-6 km2; the 4-12 prior piled up at its upper bound
    surv_out       = c("unif", 0.4, 0.97),
    surv_disp_rho  = c("unif", 4.5, 7),    # rates file 5.8; the 0-10 prior stayed at the prior
    prey_threshold = c("unif", prey_range[1], prey_range[2]),   # hindcast: a plateau (~10-30 females) needs about 14-17
    prey_months    = c("unif", 0.51, 12.49),
    start_mult     = c("unif", 1.5, 4)
  )
  if (map_type == "ndvi") p <- c(p, list(
    t_cover = c("unif", 0.3, 0.85),
    sd_hi   = c("unif", 0.05, 0.3),
    t_h     = c("unif", 0.3, 4)))
  p
}

# Named parameter list from the vector EasyABC passes (without its seed)
lynx_param_list <- function(values, map_type){
  stats::setNames(as.list(values), names(lynx_priors(map_type)))
}

# Everything a simulation needs, as plain vectors and matrices (no SpatRasters), so it can be
# sent to parallel workers.
#   template          large template (path or SpatRaster)
#   rabbit_density    monthly rabbit density stack (time = first of the month)
#   map_type          "corine": habitat_static (0/1/2 map, e.g. lynx_habitat_corine_500_large.tif)
#                     "ndvi": ndvi_mean, ndvi_sd (monthly stacks), height, barrier, roads
#   park              0/1 park map (survival inside / outside, lynx_pops_map())
#   presence_dir, counts_file   observations (read_lynx_presence(), read_lynx_counts())
#   start_year        first simulated year; the start population is put in the presence area of
#                     that year (or the latest year before it)
#   n_start           number of start lynx; default 2 x the observed breeding females of start_year
#   n_sites           number of start sites (clusters of the presence area), see lynx_start_sites()
prepare_lynx_inputs <- function(template, rabbit_density, map_type = c("corine", "ndvi"),
                                habitat_static = NULL, ndvi_mean = NULL, ndvi_sd = NULL, height = NULL,
                                barrier = NULL, roads = NULL, park, presence_dir = NULL, counts_file = NULL,
                                start_year, end_year, n_start = NULL, n_sites = 3, ndvi_climatology = FALSE){
  map_type <- match.arg(map_type)
  rd <- function(x) if (is.character(x)) terra::rast(x) else x
  template <- rd(template)[[1]]
  years <- start_year:end_year
  land <- !is.na(terra::values(template, mat = FALSE))
  on_tmpl <- function(r) terra::values(terra::mask(terra::resample(rd(r), template, method = "near"), template), mat = FALSE)

  # a monthly stack (path or SpatRaster), or read_rabbit_maps() (rabbit_maps.R): mode "mean" one
  # density for all simulations, mode "runs" one rabbit run per simulation (lynx_simulate())
  if (inherits(rabbit_density, "rabbit_maps")) {
    dens <- rabbit_density$density
    dens_months <- rabbit_density$months
  } else {
    dens <- rd(rabbit_density)
    dens_months <- as.Date(terra::time(dens))
    dens <- terra::values(dens)
    dens[!land, ] <- NA
  }

  inp <- list(map_type = map_type, years = years, ncol = terra::ncol(template), nrow = terra::nrow(template),
              land = land, template = terra::wrap(template),
              density = dens, density_months = dens_months,
              pops = on_tmpl(lynx_pops_map(template, rd(park))),
              breeding = as.numeric(land))
  if (inherits(rabbit_density, "rabbit_maps") && rabbit_density$mode == "runs") {
    inp$density_runs <- rabbit_density$runs        # land cells x months, one matrix per run
    inp$density_run_names <- rabbit_density$run_names
  }
  inp$breeding[!land] <- NA

  if (map_type == "corine") {
    inp$habitat <- on_tmpl(habitat_static)
  } else {
    nm <- rd(ndvi_mean); ns <- rd(ndvi_sd)
    # lynx years with all 12 NDVI months; with ndvi_climatology the others get the per-cell
    # median of those years (forecast years beyond the NDVI data), else they are an error
    full <- years[vapply(years, function(y) !is.null(select_year(nm, y)) && terra::nlyr(select_year(nm, y)) == 12 &&
                                                !is.null(select_year(ns, y)) && terra::nlyr(select_year(ns, y)) == 12, TRUE)]
    missing <- setdiff(years, full)
    if (length(missing) && !ndvi_climatology) stop("No full year of NDVI for the lynx year(s) ", paste(missing, collapse = ", "),
                                                   " (June of the year before to May)", call. = FALSE)
    if (length(full) == 0) stop("No lynx year with a full year of NDVI", call. = FALSE)
    yearly <- function(r, fun) vapply(full, function(y) terra::values(terra::app(select_year(r, y), fun, na.rm = TRUE), mat = FALSE), numeric(length(land)))
    nmax <- yearly(nm, "max"); smean <- yearly(ns, "mean")
    fill <- function(m) {
      out <- matrix(NA_real_, nrow(m), length(years), dimnames = list(NULL, years))
      out[, as.character(full)] <- m
      if (length(missing)) out[, as.character(missing)] <- apply(m, 1, stats::median, na.rm = TRUE)
      out
    }
    inp$ndvi_max <- fill(nmax); inp$ndvi_sdm <- fill(smean)
    inp$ndvi_climatology_years <- missing
    inp$height  <- on_tmpl(height)
    inp$barrier <- on_tmpl(barrier)
    inp$roads   <- on_tmpl(roads)
  }

  if (!is.null(presence_dir)) inp$presence <- read_lynx_presence(presence_dir, template)
  if (!is.null(counts_file))  inp$counts   <- read_lynx_counts(counts_file)

  # start sites: n_sites clusters of the presence area of the start year (or the latest year
  # before it), on CORINE-or-NDVI-independent cells: towns and roads excluded
  cells <- if (!is.null(inp$presence)) {
    py <- as.integer(colnames(inp$presence))
    sy <- max(py[py <= start_year], min(py))
    which(inp$presence[, as.character(sy)] %in% 1)
  } else which(land & inp$pops == 4)
  if (map_type == "corine") {
    ok <- cells[inp$habitat[cells] %in% 2]
  } else {
    ok <- cells[!(inp$barrier[cells] %in% 1) & !(inp$roads[cells] %in% 1)]
  }
  if (length(ok) > 0) cells <- ok
  inp$start_sites <- lynx_start_sites(cells, template, n_sites)
  if (is.null(n_start)) {
    n_start <- if (!is.null(inp$counts) && start_year %in% inp$counts$year)
      2 * inp$counts$breeding_females[inp$counts$year == start_year] else 30
  }
  inp$n_start <- round(n_start)
  # start_mult (a calibrated parameter) replaces n_start for parameter sets that have it
  inp$start_females <- if (!is.null(inp$counts) && start_year %in% inp$counts$year)
    inp$counts$breeding_females[inp$counts$year == start_year] else inp$n_start / 2
  inp
}

# Cells of the outer `width` rows and columns of the map
map_frame_cells <- function(ncol, nrow, width = 2){
  rc <- expand.grid(col = seq_len(ncol), row = seq_len(nrow))
  which(rc$col <= width | rc$col > ncol - width | rc$row <= width | rc$row > nrow - width)
}

# Habitat values (0/1/2) of every year: matrix cells x years.
# The outer two rows and columns are barrier: the IBM looks at the neighbours of a lynx's cell and
# of territory cells without checking the map edge (the PVA maps had sea around them), so a lynx
# on the last row or column read outside the map and crashed the program.
lynx_habitat_values <- function(inp, p){
  if (inp$map_type == "corine") {
    h <- matrix(inp$habitat, length(inp$habitat), length(inp$years), dimnames = list(NULL, inp$years))
  } else {
    cover <- (inp$ndvi_max >= p$t_cover & inp$ndvi_sdm <= p$sd_hi) | (inp$height >= p$t_h)
    cover[is.na(cover)] <- FALSE
    h <- 1 + cover * (inp$roads %in% 0)          # 2 = cover off roads, 1 = matrix
    h[inp$barrier %in% 1, ] <- 0
    h[!inp$land, ] <- NA
  }
  h[map_frame_cells(inp$ncol, inp$nrow), ] <- 0
  h
}

# Runs one simulation. Returns list(dir, territories, occupancy, habitat, prey) — dir is deleted
# unless keep_dir = TRUE. `p` is a named list (lynx_param_list()).
lynx_simulate <- function(p, inp, exe, base_rates, seed, work_dir = tempfile("lynx_"),
                          keep_dir = FALSE, write_maps = TRUE){
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  if (!keep_dir) on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  set.seed(seed)

  hab_dir <- file.path(work_dir, "habitat"); prey_dir <- file.path(work_dir, "prey")
  dir.create(hab_dir); dir.create(prey_dir)
  hab  <- lynx_habitat_values(inp, p)
  # with the maps of several rabbit runs: one run, drawn with the seed (the seed is then set
  # again, so the rest of the simulation is as with a single rabbit density)
  rabbit_run <- NULL
  if (!is.null(inp$density_runs)) {
    rabbit_run <- sample.int(length(inp$density_runs), 1)
    set.seed(seed)
    prey <- matrix(NA_real_, length(inp$land), length(inp$years), dimnames = list(NULL, inp$years))
    prey[inp$land, ] <- lynx_prey_matrix(inp$density_runs[[rabbit_run]], inp$years, p$prey_threshold,
                                         round(p$prey_months), inp$density_months)
  } else {
    prey <- lynx_prey_matrix(inp$density, inp$years, p$prey_threshold, round(p$prey_months), inp$density_months)
  }
  for (y in as.character(inp$years)) {
    write_lynx_values(hab[, y], inp$ncol, inp$nrow, file.path(hab_dir, sprintf("Lynx_HabitatMap_%s.txt", y)))
    write_lynx_values(prey[, y], inp$ncol, inp$nrow, file.path(prey_dir, sprintf("Lynx_PreyMap_%s.txt", y)))
  }
  breeding_file <- file.path(work_dir, "Lynx_BreedingMap.txt")
  pops_file     <- file.path(work_dir, "Lynx_pops.txt")
  write_lynx_values(inp$breeding, inp$ncol, inp$nrow, breeding_file)
  write_lynx_values(inp$pops, inp$ncol, inp$nrow, pops_file)

  template <- terra::unwrap(inp$template)
  start_file <- file.path(work_dir, "start_pop.txt")
  n_start <- if (!is.null(p$start_mult)) round(p$start_mult * inp$start_females) else inp$n_start
  write_lynx_start_pop(n_start, inp$start_sites, template, start_file)

  rates_file <- file.path(work_dir, "rates.txt")
  write_lynx_rates(base_rates, rates_file,
                   rates = c(list(Tsize = round(p$Tsize), surv_disp_rho = p$surv_disp_rho), p["rep_prob"][!is.null(p$rep_prob)]),
                   surv_adult = c("0" = p$surv_out, if (!is.null(p$surv_in)) c("2" = p$surv_in)))

  settings_file <- file.path(work_dir, "settings.txt")
  write_lynx_settings(list(
    seed = seed, start_year = min(inp$years), end_year = max(inp$years),
    create_additional_output = 0, all_year_maps = as.integer(write_maps), track_pedigree = 0,
    lynx_demography = rates_file, mapname_lynx = 0, breeding_file = breeding_file, prey_file = 0,
    map_lynx_pops = pops_file, lynx_start_size = start_file, lynx_reintro = 0,
    habitat_folder = hab_dir, breeding_folder = 0, prey_folder = prey_dir), settings_file)

  out_dir <- file.path(work_dir, "out")
  run_lynx_model(exe, settings_file, out_dir)
  terr <- read_lynx_territories(out_dir, years = inp$years)
  occ  <- if (write_maps) read_lynx_occupancy(out_dir, inp$years, template) else NULL
  list(dir = if (keep_dir) out_dir, territories = terr, occupancy = occ, habitat = hab, prey = prey,
       territories_by_pop = read_lynx_territories(out_dir, by_pop = TRUE),
       rabbit_run = if (!is.null(rabbit_run)) inp$density_run_names[rabbit_run])
}
