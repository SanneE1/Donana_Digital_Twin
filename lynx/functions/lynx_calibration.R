# ----------------------------------------------------------------------------
# ABC calibration of the lynx IBM and its evaluation (EasyABC, ABC-Lenormand as for the rabbit)
#
# Summary statistics (lynx_fit_stats(), target 0 for both), over fit_years:
#   d_count     RMSE of log(1 + territorial females) against log(1 + observed breeding females);
#               with late_from split in d_count_early and d_count_late (years from late_from on)
#   d_presence  1 - mean Dice overlap of simulated territories and observed presence (coarse cells);
#               left out of the ABC with use_presence = FALSE (then only a check in the evaluation)
# ----------------------------------------------------------------------------

# Model function for EasyABC: params[1] is the seed EasyABC passes (also the IBM's seed), then the
# parameters in the order of lynx_priors(map_type). Self-contained, so it can run on PSOCK workers:
# function_dir is sourced there when the functions are missing.
make_lynx_model_fn <- function(inp, exe, base_rates, fit_years, fact, function_dir, late_from = NULL,
                               use_presence = TRUE, n_reps = 1){
  force(inp); force(exe); force(base_rates); force(fit_years); force(fact); force(function_dir); force(late_from)
  force(use_presence); force(n_reps)
  n_count <- if (is.null(late_from)) 1 else 2
  fn <- function(params){
    if (!exists("lynx_simulate", mode = "function")) {
      for (f in list.files(function_dir, "[.]R$", full.names = TRUE)) source(f)
    }
    p <- lynx_param_list(params[-1], inp$map_type)
    # n_reps runs with different seeds, statistics averaged: the outcome of one parameter set
    # varies a lot between seeds (dies out or grows), and with one run the ABC accepts sets
    # that fit by chance
    template <- terra::unwrap(inp$template)
    reps <- lapply(seq_len(n_reps), function(k){
      seed <- (params[1] + (k - 1) * 7919) %% 2147483647
      tryCatch({
        s <- lynx_simulate(p, inp, exe, base_rates, seed = seed)
        st <- lynx_fit_stats(s$territories, s$occupancy, inp$counts, inp$presence, template,
                             fit_years, fact, late_from)$stats
        as.numeric(if (use_presence) st else st[names(st) != "d_presence"])
      }, error = function(e) {
        message("simulation failed: ", conditionMessage(e))
        c(rep(10, n_count), if (use_presence) 1)
      })
    })
    colMeans(do.call(rbind, reps))
  }
  environment(fn) <- list2env(list(inp = inp, exe = exe, base_rates = base_rates, fit_years = fit_years,
                                   fact = fact, function_dir = function_dir, late_from = late_from,
                                   use_presence = use_presence, n_count = n_count, n_reps = n_reps), parent = globalenv())
  fn
}

lynx_calibrate <- function(inp, exe, base_rates, fit_years, fact = 4, nb_simul = 100, p_acc_min = 0.2,
                           n_cluster = 1, function_dir, late_from = NULL, use_presence = TRUE,
                           n_reps = 1){
  priors <- lynx_priors(inp$map_type, inp$prey_range)
  min_simul <- 4 * length(priors)
  if (nb_simul < min_simul) {
    stop(sprintf("nb_simul = %d is too small for %d parameters: ABC-Lenormand needs nb_simul >= %d",
                 nb_simul, length(priors), min_simul), call. = FALSE)
  }
  if (!is.null(late_from) && !(any(fit_years < late_from) && any(fit_years >= late_from)))
    stop("late_start_year must split the fit years into two non-empty periods", call. = FALSE)
  model_fn <- make_lynx_model_fn(inp, exe, base_rates, fit_years, fact, function_dir, late_from, use_presence, n_reps)
  n_stats <- (if (is.null(late_from)) 1 else 2) + use_presence
  abc <- EasyABC::ABC_sequential(method = "Lenormand", model = model_fn, prior = unname(priors),
                                 inside_prior = TRUE, p_acc_min = p_acc_min,
                                 summary_stat_target = rep(0, n_stats), nb_simul = as.integer(nb_simul),
                                 use_seed = TRUE, verbose = FALSE, n_cluster = as.integer(n_cluster))
  post <- as.data.frame(abc$param)
  names(post) <- names(priors)
  post$weights <- abc$weights
  list(abc = abc, posterior = post)
}

# Posterior predictive simulations: n_draws parameter sets drawn by weight, each with its own
# seed, over all years of `inp`. Returns per-year territory counts, per-year occupancy frequency
# (share of runs with a territory, per cell) and the per-run statistics.
lynx_posterior_predictive <- function(posterior, inp, exe, base_rates, n_draws = 50, seed = 1,
                                      fit_years = NULL, fact = 4, late_from = NULL){
  set.seed(seed)
  idx <- sample(nrow(posterior), n_draws, replace = TRUE, prob = posterior$weights)
  pn <- setdiff(names(posterior), "weights")
  template <- terra::unwrap(inp$template)
  occ_sum <- NULL
  counts <- list(); stats <- list(); dice <- list()
  for (i in seq_len(n_draws)) {
    p <- as.list(posterior[idx[i], pn])
    s <- tryCatch(lynx_simulate(p, inp, exe, base_rates, seed = seed * 100000 + i),
                  error = function(e) { message("run ", i, " failed: ", conditionMessage(e)); NULL })
    if (is.null(s)) next
    counts[[length(counts) + 1]] <- data.frame(run = i, draw = idx[i], rabbit_run = s$rabbit_run, s$territories)
    occ_sum <- if (is.null(occ_sum)) s$occupancy else occ_sum + s$occupancy
    if (!is.null(fit_years)) {
      fs <- lynx_fit_stats(s$territories, s$occupancy, inp$counts, inp$presence, template, fit_years, fact, late_from)
      stats[[length(stats) + 1]] <- data.frame(run = i, t(fs$stats))
    }
    all_fs <- lynx_fit_stats(s$territories, s$occupancy, inp$counts, inp$presence, template, inp$years, fact)
    if (length(all_fs$dice)) dice[[length(dice) + 1]] <- data.frame(run = i, year = as.integer(names(all_fs$dice)), dice = all_fs$dice)
  }
  n_ok <- length(counts)
  list(counts = do.call(rbind, counts), stats = do.call(rbind, stats), dice = do.call(rbind, dice),
       occupancy = if (n_ok) occ_sum / n_ok, n_runs = n_ok, draws = idx)
}

# Evaluation of a posterior predictive run against the observations, split in calibration years
# and held-out years. Returns a table of metrics (one row per period and statistic) and the
# per-year table used for the plots.
lynx_evaluate <- function(pp, inp, fit_years){
  cnt <- pp$counts
  yr <- sort(unique(cnt$year))
  q <- do.call(rbind, lapply(yr, function(y){
    v <- cnt$females_territorial[cnt$year == y]
    data.frame(year = y, median = stats::median(v), lo = stats::quantile(v, 0.05), hi = stats::quantile(v, 0.95),
               p_extinct = mean(v == 0))
  }))
  q$obs <- inp$counts$breeding_females[match(q$year, inp$counts$year)]
  q$period <- ifelse(q$year %in% fit_years, "calibration", ifelse(q$year < min(fit_years), "spin-up", "held out"))
  dq <- if (!is.null(pp$dice)) stats::aggregate(dice ~ year, pp$dice, stats::median) else NULL
  q$dice_median <- if (!is.null(dq)) dq$dice[match(q$year, dq$year)] else NA

  metric <- function(period){
    d <- q[q$period == period & !is.na(q$obs), ]
    dd <- q[q$period == period & !is.na(q$dice_median), ]
    data.frame(period = period,
               n_count_years = nrow(d),
               rmse_log_count = if (nrow(d)) sqrt(mean((log1p(d$median) - log1p(d$obs))^2)) else NA,
               bias_count = if (nrow(d)) mean(d$median - d$obs) else NA,
               coverage90_count = if (nrow(d)) mean(d$obs >= d$lo & d$obs <= d$hi) else NA,
               n_presence_years = nrow(dd),
               mean_dice = if (nrow(dd)) mean(dd$dice_median) else NA,
               p_extinct_last = utils::tail(q$p_extinct, 1))
  }
  list(metrics = rbind(metric("calibration"), metric("held out")), by_year = q)
}

lynx_eval_plots <- function(ev, title = ""){
  q <- ev$by_year
  p1 <- ggplot2::ggplot(q, ggplot2::aes(year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#c8a041", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = median), colour = "#8b6b13", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = obs, shape = period), na.rm = TRUE) +
    ggplot2::scale_shape_manual(values = c(calibration = 16, "held out" = 1, "spin-up" = 4)) +
    ggplot2::labs(x = NULL, y = "Territorial / breeding females", shape = "Observed",
                  title = title, subtitle = "Line and band: posterior predictive median and 90% interval") +
    ggplot2::theme_minimal(base_size = 12)
  p2 <- ggplot2::ggplot(q[!is.na(q$dice_median), ], ggplot2::aes(year, dice_median, colour = period)) +
    ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::ylim(0, 1) +
    ggplot2::labs(x = NULL, y = "Dice overlap with observed presence (median)", colour = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12)
  list(counts = p1, presence = p2)
}
