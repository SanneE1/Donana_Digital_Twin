# ----------------------------------------------------------------------------
# Lynx calibration from file paths to results (used by the donana-lynx-calibration component and
# by lynx/Workflow_lynx_calibration.R). Nothing is written here: main.r saves the results.
# ----------------------------------------------------------------------------

lynx_calibration_run <- function(
    template, rabbit_density, rabbit_maps = c("mean", "runs"), habitat_map = c("corine", "ndvi"),
    habitat_corine = NULL, ndvi_mean = NULL, ndvi_sd = NULL, vegetation_height = NULL,
    barrier = NULL, roads = NULL, park,
    presence_dir, counts_file, rates_file, model_exe,
    start_year = 2006, end_year = 2024, fit_years = 2007:2018, late_start_year = 2012,
    use_presence = TRUE, n_replicates = 1,
    n_start = NULL, n_start_sites = 3, presence_cell_factor = 4, prey_range = NULL,
    nb_simul = 100, p_acc_min = 0.2, n_cluster = 1, n_sim_samples = 50,
    function_dir){

  habitat_map <- match.arg(habitat_map)
  rabbit_maps <- match.arg(rabbit_maps)
  if (habitat_map == "corine" && is.null(habitat_corine)) stop("habitat_map 'corine' needs the CORINE habitat map")
  if (habitat_map == "ndvi" && (is.null(ndvi_mean) || is.null(ndvi_sd) || is.null(vegetation_height) ||
                                is.null(barrier) || is.null(roads))) {
    stop("habitat_map 'ndvi' needs NDVI mean and SD, vegetation height, barrier and roads maps")
  }
  if (!all(fit_years %in% start_year:end_year)) stop("fit_years must lie within start_year:end_year")

  cat("\nPreparing inputs (", habitat_map, " habitat map, ", start_year, "-", end_year, ")\n", sep = "")
  # rabbit maps: a file or a directory with the maps per run and/or their mean (rabbit_maps.R);
  # only the months of the lynx years are read
  if (is.character(rabbit_density)) {
    rabbit_density <- read_rabbit_maps(rabbit_density, template, rabbit_maps,
                                       months = lynx_years_months(start_year, end_year))
    cat(sprintf("rabbit maps: %s (%s)\n", rabbit_density$source, rabbit_maps))
  } else if (rabbit_maps == "runs" && !inherits(rabbit_density, "rabbit_maps")) {
    stop("rabbit_maps = 'runs' needs the path of the rabbit maps per run")
  }
  inp <- prepare_lynx_inputs(template, rabbit_density, habitat_map,
                             habitat_static = habitat_corine, ndvi_mean = ndvi_mean, ndvi_sd = ndvi_sd,
                             height = vegetation_height, barrier = barrier, roads = roads, park = park,
                             presence_dir = presence_dir, counts_file = counts_file,
                             start_year = start_year, end_year = end_year,
                             n_start = n_start, n_sites = n_start_sites)
  cat(sprintf("start population: %d lynx at %d sites; presence years %s; count years %s\n",
              inp$n_start, length(inp$start_sites), paste(range(colnames(inp$presence)), collapse = "-"),
              paste(range(inp$counts$year), collapse = "-")))

  # prior of prey_threshold on the scale of the rabbit maps, unless given
  if (is.null(prey_range) || any(is.na(prey_range))) {
    auto <- prey_threshold_range(inp)
    if (is.null(prey_range)) prey_range <- c(NA, NA)
    prey_range[is.na(prey_range)] <- auto[is.na(prey_range)]
  }
  if (prey_range[1] <= 0 || prey_range[2] <= prey_range[1]) stop("prey_threshold prior: need 0 < min < max")
  inp$prey_range <- prey_range
  cat(sprintf("prey_threshold prior: %.3g-%.3g rabbits per cell\n", prey_range[1], prey_range[2]))

  model_exe  <- normalizePath(model_exe)
  rates_file <- normalizePath(rates_file)
  # A model without the Doñana-twin outputs (lynx_territories.csv) fails here, not in every run
  check <- tryCatch(lynx_simulate(lapply(lynx_priors(habitat_map, inp$prey_range), function(x) mean(as.numeric(x[2:3]))),
                                  inp, model_exe, rates_file, seed = 1),
                    error = function(e) stop("Test run of the lynx model failed: ", conditionMessage(e), call. = FALSE))

  cat("\nRunning ABC\n")
  t0 <- Sys.time()
  cal <- lynx_calibrate(inp, model_exe, rates_file, fit_years, presence_cell_factor,
                        nb_simul = nb_simul, p_acc_min = p_acc_min, n_cluster = n_cluster,
                        function_dir = function_dir, late_from = late_start_year,
                        use_presence = use_presence, n_reps = n_replicates)
  cat(sprintf("ABC finished: %d simulations in %.1f min\n", cal$abc$nsim,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  cat("\nPosterior predictive evaluation\n")
  pp <- lynx_posterior_predictive(cal$posterior, inp, model_exe, rates_file, n_draws = n_sim_samples,
                                  fit_years = fit_years, fact = presence_cell_factor,
                                  late_from = late_start_year)
  ev <- lynx_evaluate(pp, inp, fit_years)
  print(ev$metrics, digits = 3)
  plots <- lynx_eval_plots(ev, sprintf("Doñana lynx, %s habitat map", habitat_map))

  # share of posterior predictive runs with a territory, per cell and year
  occ <- terra::rast(lapply(seq_len(ncol(pp$occupancy)), function(j) cells_to_map(pp$occupancy[, j], terra::unwrap(inp$template))))
  names(occ) <- colnames(pp$occupancy)

  priors <- lynx_priors(habitat_map, inp$prey_range)
  summary <- data.frame(parameter = names(priors),
                        prior_lo = sapply(priors, function(x) as.numeric(x[2])),
                        prior_hi = sapply(priors, function(x) as.numeric(x[3])),
                        post_mean = sapply(names(priors), function(n) stats::weighted.mean(cal$posterior[[n]], cal$posterior$weights)),
                        post_sd = sapply(names(priors), function(n) {
                          w <- cal$posterior$weights / sum(cal$posterior$weights); m <- sum(w * cal$posterior[[n]])
                          sqrt(sum(w * (cal$posterior[[n]] - m)^2)) }),
                        row.names = NULL)
  summary$sd_ratio <- summary$post_sd / ((summary$prior_hi - summary$prior_lo) / sqrt(12))

  list(posterior = cal$posterior, abc = cal$abc, metrics = cbind(habitat_map = habitat_map, ev$metrics),
       by_year = ev$by_year, runs = pp$counts, run_stats = pp$stats, plots = plots, occupancy = occ,
       parameter_summary = summary)
}
