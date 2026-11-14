# ----------------------------------------------------------------------------
# Lynx forecast: an ensemble of IBM runs with the calibrated posterior
#
# The IBM cannot start from a saved population state, so every run starts in start_year (as in
# the calibration, start population in the observed presence area) and runs through the
# observed years into the forecast years: the population at the forecast start then carries the
# history of the hindcast. Rabbit density comes from one or more member files (monthly stacks on
# the template, time = first of the month) that cover start_year - 1 (June) to end_year (May): e.g.
# a rabbit hindcast followed by a rabbit forecast. Each run uses one member.
#
# Uncertainty sources (each can be switched off):
#   parameters   a posterior parameter set drawn by weight   (off: weighted posterior mean)
#   demographic  the IBM's random seed                        (off: one seed for every run)
#   rabbit       a rabbit density member                      (off: the mean of the members)
#
# The NDVI habitat map needs NDVI for every lynx year; years beyond the NDVI data use the per-cell
# median of the yearly NDVI values (climatology).
# ----------------------------------------------------------------------------

# Rabbit density members: files (each a monthly stack) -> list of matrices cells x months
read_density_members <- function(files, template){
  lapply(files, function(f){
    r <- terra::rast(f)
    t <- as.Date(terra::time(r))
    if (any(is.na(t))) stop("Rabbit density file without time per layer: ", f)
    if (!terra::compareGeom(r, template, stopOnError = FALSE)) r <- terra::resample(r, template, method = "near")
    list(values = terra::values(r), months = t)
  })
}

lynx_forecast_run <- function(
    template, density_files, habitat_map = c("corine", "ndvi"),
    habitat_corine = NULL, ndvi_mean = NULL, ndvi_sd = NULL, vegetation_height = NULL,
    barrier = NULL, roads = NULL, park, presence_dir, counts_file = NULL,
    posterior_file, rates_file, model_exe,
    start_year = 2006, forecast_start_year, end_year,
    n_runs = 100, uncertainty_sources = c("parameters", "demographic", "rabbit"),
    n_start = NULL, n_start_sites = 3, seed = 1){

  habitat_map <- match.arg(habitat_map)
  bad <- setdiff(uncertainty_sources, c("parameters", "demographic", "rabbit"))
  if (length(bad)) stop("Unknown uncertainty source(s): ", paste(bad, collapse = ", "),
                        " (use parameters, demographic, rabbit)", call. = FALSE)
  template <- terra::rast(template)[[1]]
  set.seed(seed)

  members <- read_density_members(density_files, template)
  months_needed <- seq(as.Date(sprintf("%d-06-01", start_year - 1)), as.Date(sprintf("%d-05-01", end_year)), by = "month")
  for (i in seq_along(members)) {
    miss <- setdiff(format(months_needed, "%Y-%m"), format(members[[i]]$months, "%Y-%m"))
    if (length(miss)) stop("Rabbit density member ", basename(density_files[i]), " lacks ", length(miss),
                           " month(s), e.g. ", paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)
  }
  # common months, mean of the members
  common <- Reduce(intersect, lapply(members, function(m) format(m$months, "%Y-%m")))
  members <- lapply(members, function(m){ k <- match(common, format(m$months, "%Y-%m"));
                                          list(values = m$values[, k, drop = FALSE], months = m$months[k]) })
  dens_mean <- Reduce(`+`, lapply(members, `[[`, "values")) / length(members)
  mean_rast <- terra::rast(template, nlyrs = length(common))
  terra::values(mean_rast) <- dens_mean
  terra::time(mean_rast) <- members[[1]]$months

  inp <- prepare_lynx_inputs(template, mean_rast, habitat_map,
                             habitat_static = habitat_corine, ndvi_mean = ndvi_mean, ndvi_sd = ndvi_sd,
                             height = vegetation_height, barrier = barrier, roads = roads, park = park,
                             presence_dir = presence_dir, counts_file = counts_file,
                             start_year = start_year, end_year = end_year,
                             n_start = n_start, n_sites = n_start_sites, ndvi_climatology = TRUE)

  post <- utils::read.csv(posterior_file)
  pn <- names(lynx_priors(habitat_map))
  if (!all(pn %in% names(post))) stop("The posterior lacks parameter(s) ", paste(setdiff(pn, names(post)), collapse = ", "),
                                      " of the ", habitat_map, " habitat map: was it calibrated with another map, or before ",
                                      "the parameter was added?", call. = FALSE)
  pn <- setdiff(names(post), "weights")   # incl. rep_prob / surv_in of posteriors that calibrated them
  w <- post$weights / sum(post$weights)
  runs <- data.frame(run = seq_len(n_runs),
                     draw = if ("parameters" %in% uncertainty_sources) sample(nrow(post), n_runs, TRUE, prob = w) else NA,
                     seed = if ("demographic" %in% uncertainty_sources) sample.int(1e8, n_runs) else seed,
                     member = if ("rabbit" %in% uncertainty_sources) rep_len(sample(length(members)), n_runs) else NA)
  mean_par <- as.list(colSums(post[, pn] * w))

  exe <- normalizePath(model_exe); rates_file <- normalizePath(rates_file)
  counts <- list(); occ_sum <- 0; n_ok <- 0
  for (i in seq_len(n_runs)) {
    p <- if (is.na(runs$draw[i])) mean_par else as.list(post[runs$draw[i], pn])
    inp_i <- inp
    if (!is.na(runs$member[i])) inp_i$density <- members[[runs$member[i]]]$values
    s <- tryCatch(lynx_simulate(p, inp_i, exe, rates_file, seed = runs$seed[i]),
                  error = function(e) { message("run ", i, " failed: ", conditionMessage(e)); NULL })
    if (is.null(s)) next
    counts[[length(counts) + 1]] <- data.frame(run = i, s$territories)
    occ_sum <- occ_sum + s$occupancy; n_ok <- n_ok + 1
    if (i %% 10 == 0) cat("run", i, "of", n_runs, "\n")
  }
  if (n_ok == 0) stop("All forecast runs failed")
  counts <- do.call(rbind, counts)

  yr <- sort(unique(counts$year))
  # latest observed breeding females before the forecast start
  before <- if (!is.null(inp$counts)) inp$counts[inp$counts$year < forecast_start_year, ] else NULL
  last_obs <- if (!is.null(before) && nrow(before)) before$breeding_females[which.max(before$year)] else numeric(0)
  summ <- do.call(rbind, lapply(yr, function(y){
    v <- counts$females_territorial[counts$year == y]
    a <- counts$all_lynx[counts$year == y]
    data.frame(year = y, period = if (y < forecast_start_year) "hindcast" else "forecast",
               females_median = stats::median(v), females_lo90 = stats::quantile(v, 0.05), females_hi90 = stats::quantile(v, 0.95),
               females_lo50 = stats::quantile(v, 0.25), females_hi50 = stats::quantile(v, 0.75),
               all_lynx_median = stats::median(a),
               p_extinct = mean(v == 0),
               p_below_last_observed = if (length(last_obs)) mean(v < last_obs) else NA)
  }))
  summ$observed <- if (!is.null(inp$counts)) inp$counts$breeding_females[match(summ$year, inp$counts$year)] else NA

  occ <- terra::rast(lapply(seq_len(ncol(occ_sum)), function(j) cells_to_map(occ_sum[, j] / n_ok, template)))
  names(occ) <- colnames(occ_sum)

  plot <- ggplot2::ggplot(summ, ggplot2::aes(year)) +
    ggplot2::annotate("rect", xmin = forecast_start_year - 0.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey90") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = females_lo90, ymax = females_hi90), fill = "#c8a041", alpha = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = females_lo50, ymax = females_hi50), fill = "#c8a041", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = females_median), colour = "#8b6b13", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = observed), na.rm = TRUE) +
    ggplot2::labs(x = NULL, y = "Territorial / breeding females",
                  title = sprintf("Doñana lynx forecast (%s habitat map), %d runs", habitat_map, n_ok),
                  subtitle = sprintf("Bands: 50%% and 90%% of runs; points: observed; grey: forecast from %d. Uncertainty: %s",
                                     forecast_start_year, paste(uncertainty_sources, collapse = ", "))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"))

  list(summary = summ, runs = counts, run_table = runs, occupancy = occ, plot = plot,
       metadata = list(habitat_map = habitat_map, start_year = start_year, forecast_start_year = forecast_start_year,
                       end_year = end_year, n_runs = n_runs, n_ok = n_ok, uncertainty_sources = uncertainty_sources,
                       density_members = basename(density_files), seed = seed))
}
