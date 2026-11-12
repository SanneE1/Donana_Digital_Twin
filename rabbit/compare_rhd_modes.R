#  Which RHD mode fits the observed rabbit KAI best: none, schedule or immunity?
#
#  The protocol is in rabbit/Pascal_program/model_settings.md ("Deciding the
#  default RHD mode"). In short, for every mode with the same data, priors,
#  ABC settings and seeds:
#    1. calibrate on the two shared KAI statistics only (median r and median r
#       of the changes over the transects), for 2017-2024 and for 2017-2023
#    2. posterior predictive fit: n_draws parameter sets x seeds, in-sample
#    3. held-out year: the 2017-2023 posteriors predict the 2024 censuses
#    4. demographic noise: n_noise seeds at the weighted posterior mean
#  and compare the modes with bootstrap intervals of the differences.
#  Immunity becomes the default only if it beats schedule in-sample (interval
#  excludes 0) and is not worse on the held-out year.
#
#  It uses the calibration component's code (rabbit_calibration(), model_fn()),
#  so the comparison runs exactly the code that calibrates. Needs the rabbit
#  model build with named arguments (branch named-model-settings).
#
#  Run from the Donana_Digital_Twin folder:  Rscript rabbit/compare_rhd_modes.R [quick|small]
#  The full run takes many hours (6 calibrations with nb_simul = 300, and about
#  3 x (n_draws x seeds x 2 + n_noise) simulations); "small" is a reduced run
#  (about a quarter of the work, wider intervals) that fits in one night on a
#  laptop; "quick" is a smoke test only.
#  Results: out_dir/ (posteriors, simulation scores, summary.csv, decision.txt).

args  <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
small <- !quick && "small" %in% args
level <- if (quick) "quick" else if (small) "small" else "full"

# ---- Settings -------------------------------------------------------------------------------

components <- normalizePath(file.path("..", "Components"), winslash = "/")
comp_code  <- file.path(components, "donana-rabbit-calibration", "code")
data       <- file.path(components, "data")
model_exe  <- normalizePath(file.path("rabbit", "Pascal_program", "programs", "Rabbit_x86_64-win64.exe"))
rates      <- normalizePath(file.path("data", "rabbit", "Rabbit_demography_values.txt"))
out_dir    <- file.path("rabbit", "results",
                        if (level == "full") "rhd_mode_comparison" else paste0("rhd_mode_comparison_", level))

inputs_args <- list(
  template            = file.path(data, "template_raster_500.tif"),
  obs_path            = file.path(data, "rabbit_census", "KAI_Rabbit_Night_2024_v1.csv"),
  temperature_files   = file.path(data, "adapted", "temperature"),
  is_Kelvin           = TRUE,
  precipitation_files = file.path(data, "adapted", "precipitation"),
  ndvi_max_files      = file.path(data, "adapted", "ndvi_filled"),
  ndvi_sd_files       = file.path(data, "adapted", "ndvi_sd_filled"),
  flood_files         = file.path(data, "adapted", "flooding"),
  height_files        = file.path(data, "height.tif"),
  fill_gaps           = TRUE)
transects <- normalizePath(file.path(data, "rabbit_census", "Transect_oryctolagus.kml"))

modes        <- c("none", "schedule", "immunity")
start_year   <- 2017
end_year     <- 2024     # in-sample period end; the held-out calibration ends a year earlier
#                        quick  small   full
nb_simul     <- switch(level, quick = 68,  small = 150, full = 300)
p_acc_min    <- switch(level, quick = 0.6, small = 0.2, full = 0.05)
n_draws      <- switch(level, quick = 4,   small = 100, full = 200)
seeds        <- switch(level, quick = 1,   small = 1:2, full = 1:3)
n_noise      <- switch(level, quick = 3,   small = 20,  full = 50)
n_cluster    <- switch(level, quick = 2,   small = 3,   full = 4)
# Workers can override the default:  Rscript rabbit/compare_rhd_modes.R small cores=2
cores_arg <- grep("^cores=", args, value = TRUE)
if (length(cores_arg)) n_cluster <- as.integer(sub("^cores=", "", cores_arg[1]))
n_boot       <- 2000

PARAMS13 <- c("K_max", "ndvi_opt", "ndvi_sigma", "ndvi_sd_opt", "ndvi_sd_sigma", "height_opt",
              "height_sigma", "dens_opt", "R_lambda", "R_sigma", "ndvi_w", "ndvi_sd_w", "height_w")

# ---- Setup ----------------------------------------------------------------------------------

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(dplyr); library(tidyr); library(lubridate)
  library(readr); library(sf); library(EasyABC); library(ggplot2)
})
for (f in list.files(comp_code, pattern = "\\.R$", full.names = TRUE)) source(f)
# rabbit_calibration() reads these globals (set by the component's main.r):
# the parallel workers source the code from function_dir and start in proj_root
function_dir <- normalizePath(comp_code, winslash = "/")
proj_root    <- normalizePath(getwd(), winslash = "/")
check_rabbit_model(model_exe)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Calibrations --------------------------------------------------------------------------

calibrate <- function(mode, last_year) {
  file <- file.path(out_dir, sprintf("posterior_%s_%d-%d.csv", mode, start_year, last_year))
  if (file.exists(file)) {
    message("Using ", file)
    return(read.csv(file))
  }
  set.seed(1)
  res <- do.call(rabbit_calibration, c(inputs_args, list(
    start_year = start_year, end_year = last_year, obs_transects = transects,
    rabbit_demographic_rates = rates, demograph_mod = model_exe,
    N_sim_samples = 2, n_cluster = n_cluster, nb_simul = nb_simul, p_acc_min = p_acc_min,
    rhd_mode = mode, rhd_statistic = FALSE)))
  write.csv(res$posterior, file, row.names = FALSE)
  res$posterior
}

posteriors <- list()
for (mode in modes) {
  posteriors[[paste(mode, "full")]]    <- calibrate(mode, end_year)
  posteriors[[paste(mode, "heldout")]] <- calibrate(mode, end_year - 1)
}

# ---- Simulations with the same inputs -------------------------------------------------------

clim_dir <- tempfile("clim_compare")
dir.create(clim_dir)
inputs <- do.call(prepare_calibration_inputs, c(inputs_args, list(
  start_year = start_year, end_year = end_year, out_dir_clim = clim_dir)))

# One simulation: the parameter vector model_fn() takes (seed, 13 parameters,
# RHD parameters of the mode) -> simulated and observed KAI per census
sim_state <- list(
  Cov_list = lapply(inputs$Cov_list, terra::wrap), template = terra::wrap(inputs$DT_template),
  obs_df = inputs$obs_df, clim = normalizePath(clim_dir), exe = model_exe, rates = rates,
  transects = transects, start_year = start_year, end_year = end_year, comp_code = comp_code)

simulate_one <- function(job, st) {
  setup <- rhd_calibration_setup(job$mode)
  out <- model_fn(params = job$params, start_yr = st$start_year, end_yr = st$end_year,
                  Cov_list = lapply(st$Cov_list, terra::unwrap),
                  DT_template = terra::unwrap(st$template),
                  demograph_mod = st$exe, rabbit_demographic_rates = st$rates,
                  obs_transects = st$transects, obs_df = st$obs_df, out_dir_clim = st$clim,
                  rhd_settings = setup$settings, rhd_param_names = names(setup$priors),
                  return_sim = TRUE)
  if (!is.data.frame(out)) return(NULL)   # failed run
  out$job <- job$id
  out
}

simulate_jobs <- function(jobs) {
  cl <- parallel::makeCluster(n_cluster)
  on.exit(parallel::stopCluster(cl))
  init <- function(st) {
    suppressPackageStartupMessages({
      library(terra); library(dplyr); library(tidyr); library(lubridate); library(readr)
    })
    for (f in list.files(st$comp_code, pattern = "\\.R$", full.names = TRUE)) source(f)
    assign(".st", st, envir = globalenv())
    NULL
  }
  run <- function(job) simulate_one(job, get(".st", envir = globalenv()))
  environment(init) <- globalenv(); environment(run) <- globalenv()
  environment(simulate_one) <- globalenv()
  parallel::clusterExport(cl, "simulate_one", envir = environment())
  parallel::clusterCall(cl, init, sim_state)
  bind_rows(parallel::parLapplyLB(cl, jobs, run))
}

param_vector <- function(post_row, mode, seed) {
  rhd_names <- names(rhd_calibration_setup(mode)$priors)
  c(seed, unlist(post_row[c(PARAMS13, rhd_names)]))
}

weighted_mean_row <- function(post) {
  w <- post$weights / sum(post$weights)
  as.data.frame(t(colSums(post[setdiff(names(post), "weights")] * w)))
}

# ---- Scores ---------------------------------------------------------------------------------

# In-sample statistics of one simulation, as in model_fn()
insample_scores <- function(m, last_year) {
  m %>%
    filter(year(Fecha) <= last_year) %>%
    group_by(transect) %>%
    arrange(Fecha, .by_group = TRUE) %>%
    summarise(r = suppressWarnings(cor(KAI_obs, KAI_sim, use = "complete.obs")),
              r_delta = suppressWarnings(cor(diff(KAI_obs), diff(KAI_sim), use = "complete.obs")),
              rmse_scaled = sqrt(mean((KAI_obs - KAI_sim * mean(KAI_obs, na.rm = TRUE) /
                                         mean(KAI_sim, na.rm = TRUE))^2, na.rm = TRUE)) /
                            sd(KAI_obs, na.rm = TRUE),
              .groups = "drop") %>%
    summarise(r = median(r, na.rm = TRUE), r_delta = median(r_delta, na.rm = TRUE),
              rmse_scaled = mean(rmse_scaled[is.finite(rmse_scaled)]))
}

# Held-out year: simulated KAI rescaled per transect with the ratio in the
# calibration years; scaled error, sign of the change from the last
# calibration year, and correlation across transects in the held-out year
heldout_scores <- function(m, last_year) {
  rate <- m %>% filter(year(Fecha) <= last_year) %>% group_by(transect) %>%
    summarise(rate = mean(KAI_obs, na.rm = TRUE) / mean(KAI_sim, na.rm = TRUE), .groups = "drop")
  y <- m %>% left_join(rate, by = "transect") %>%
    mutate(sim = ifelse(is.finite(rate), KAI_sim * rate, NA), year = year(Fecha)) %>%
    filter(year %in% c(last_year, last_year + 1)) %>%
    group_by(transect, year) %>%
    summarise(obs = mean(KAI_obs, na.rm = TRUE), sim = mean(sim, na.rm = TRUE), .groups = "drop")
  sd_obs <- m %>% filter(year(Fecha) <= last_year) %>% group_by(transect) %>%
    summarise(sd = sd(KAI_obs, na.rm = TRUE), .groups = "drop")
  wide <- y %>% pivot_wider(names_from = year, values_from = c(obs, sim))
  names(wide) <- sub(paste0("_", last_year + 1, "$"), "_new", sub(paste0("_", last_year, "$"), "_old", names(wide)))
  wide <- wide %>% left_join(sd_obs, by = "transect")
  tibble(
    error_scaled = mean(abs(wide$sim_new - wide$obs_new) / wide$sd, na.rm = TRUE),
    sign_right   = mean(sign(wide$sim_new - wide$sim_old) == sign(wide$obs_new - wide$obs_old), na.rm = TRUE),
    r_across     = suppressWarnings(cor(wide$sim_new, wide$obs_new, use = "complete.obs")))
}

score_runs <- function(sims, jobs, fun, last_year) {
  meta <- bind_rows(lapply(jobs, function(j) tibble(job = j$id, mode = j$mode, kind = j$kind)))
  sims %>% group_by(job) %>% group_modify(~ fun(.x, last_year)) %>% ungroup() %>%
    left_join(meta, by = "job")
}

# ---- 2.-4. Jobs and simulations ---------------------------------------------------------------

make_jobs <- function(kind, post_key, mode, rows, run_seeds, start_id) {
  post <- posteriors[[post_key]]
  jobs <- list()
  for (r in rows) for (s in run_seeds) {
    row <- if (identical(r, "mean")) weighted_mean_row(post) else post[r, ]
    jobs[[length(jobs) + 1]] <- list(id = start_id + length(jobs), kind = kind, mode = mode,
                                     params = param_vector(row, mode, s))
  }
  jobs
}

set.seed(42)
jobs <- list()
for (mode in modes) {
  for (period in c("full", "heldout")) {
    key  <- paste(mode, period)
    rows <- sample(nrow(posteriors[[key]]), n_draws, replace = TRUE, prob = posteriors[[key]]$weights)
    jobs <- c(jobs, make_jobs(period, key, mode, rows, seeds, length(jobs)))
  }
  jobs <- c(jobs, make_jobs("noise", paste(mode, "full"), mode, list("mean"), 1000 + seq_len(n_noise),
                            length(jobs)))
}
for (i in seq_along(jobs)) jobs[[i]]$id <- i
message(length(jobs), " simulations")
sims <- simulate_jobs(jobs)
write_csv(sims, file.path(out_dir, "simulations.csv"))

kinds <- vapply(jobs, `[[`, character(1), "kind")
insample <- score_runs(sims %>% filter(job %in% which(kinds %in% c("full", "noise"))),
                       jobs, insample_scores, end_year)
heldout  <- score_runs(sims %>% filter(job %in% which(kinds == "heldout")),
                       jobs, heldout_scores, end_year - 1)
write_csv(insample, file.path(out_dir, "scores_insample.csv"))
write_csv(heldout, file.path(out_dir, "scores_heldout.csv"))

# ---- Comparison -------------------------------------------------------------------------------

boot_diff <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  d <- replicate(n_boot, mean(sample(a, replace = TRUE)) - mean(sample(b, replace = TRUE)))
  c(diff = mean(a) - mean(b), lo = unname(quantile(d, 0.025)), hi = unname(quantile(d, 0.975)))
}

pairs <- list(c("immunity", "schedule"), c("immunity", "none"), c("schedule", "none"))
summary_rows <- list()
for (p in pairs) {
  for (stat in c("r", "r_delta", "rmse_scaled")) {
    a <- insample %>% filter(kind == "full", mode == p[1]) %>% pull(stat)
    b <- insample %>% filter(kind == "full", mode == p[2]) %>% pull(stat)
    summary_rows[[length(summary_rows) + 1]] <- c(list(comparison = paste(p, collapse = " - "),
      set = "in-sample", statistic = stat), as.list(boot_diff(a, b)))
  }
  for (stat in c("error_scaled", "sign_right", "r_across")) {
    a <- heldout %>% filter(mode == p[1]) %>% pull(stat)
    b <- heldout %>% filter(mode == p[2]) %>% pull(stat)
    summary_rows[[length(summary_rows) + 1]] <- c(list(comparison = paste(p, collapse = " - "),
      set = paste("held-out", end_year), statistic = stat), as.list(boot_diff(a, b)))
  }
}
comparison <- bind_rows(summary_rows)
noise <- insample %>% filter(kind == "noise") %>% group_by(mode) %>%
  summarise(across(c(r, r_delta, rmse_scaled), ~ sd(.x, na.rm = TRUE), .names = "sd_{.col}"))
write_csv(comparison, file.path(out_dir, "summary.csv"))
write_csv(noise, file.path(out_dir, "demographic_noise.csv"))
print(comparison, n = Inf)
print(noise)

# Decision rule (higher r / r_delta / sign_right / r_across is better; lower errors are better)
pick <- function(which_set, which_stat) {
  comparison[comparison$comparison == "immunity - schedule" &
               startsWith(comparison$set, which_set) & comparison$statistic == which_stat, ]
}
better_in <- isTRUE(pick("in-sample", "r")$lo > 0 && pick("in-sample", "r_delta")$lo > 0)
not_worse_out <- isTRUE(pick("held-out", "error_scaled")$lo <= 0)
decision <- if (better_in && not_worse_out) "immunity" else "schedule"
msg <- sprintf(paste0("Default RHD mode: %s\n  immunity better than schedule in-sample (r and r_delta, ",
                      "interval excludes 0): %s\n  immunity not worse on the held-out year (scaled error): %s\n",
                      "Settings: nb_simul = %d, p_acc_min = %s, %d draws x %d seeds, %d noise seeds%s\n"),
               decision, better_in, not_worse_out, nb_simul, p_acc_min, n_draws, length(seeds), n_noise,
               switch(level, quick = "  (QUICK SMOKE TEST: not a valid result)",
                      small = "  (SMALL RUN: indicative; confirm with the full run)", full = ""))
cat("\n", msg)
writeLines(msg, file.path(out_dir, "decision.txt"))
unlink(clim_dir, recursive = TRUE)
