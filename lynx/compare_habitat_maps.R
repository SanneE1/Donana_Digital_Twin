#
# Compares the fit of the lynx IBM with the CORINE and the NDVI habitat map
#
# Takes the posterior of a calibration with each map (donana-lynx-calibration, posterior_samples.csv)
# and runs the same posterior predictive evaluation for both: n_draws runs (posterior draw by
# weight + seed), 2006-2024, metrics for the calibration years (2007-2018) and the held-out years
# (counts 2019-2024, presence 2020-2022), see lynx_evaluate().
#
#   Rscript lynx/compare_habitat_maps.R <corine posterior.csv> <ndvi posterior.csv> [out_dir] [n_draws]
#
# Run from the Donana_Digital_Twin project root. Writes to out_dir (default
# lynx/results/habitat_map_comparison): comparison_metrics.csv, by_year.csv, comparison.png.

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({library(terra); library(ggplot2)})
for (f in list.files(file.path("lynx", "functions"), "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Give the CORINE and NDVI posterior files")
posteriors <- c(corine = args[1], ndvi = args[2])
out_dir <- if (length(args) >= 3) args[3] else file.path("lynx", "results", "habitat_map_comparison")
n_draws <- if (length(args) >= 4) as.integer(args[4]) else 50
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ed <- file.path("data", "environmental_data")
start_year <- 2006; end_year <- 2024; fit_years <- 2007:2018
exe   <- normalizePath(file.path("lynx", "Pascal_program", "Executables", "Run_model_release.exe"))
rates <- normalizePath(file.path("data", "lynx", "Lynx_demography_values_donana.txt"))

res <- lapply(names(posteriors), function(map_type){
  inp <- prepare_lynx_inputs(file.path(ed, "template_raster_500_large.tif"),
                             file.path("data", "rabbit_hindcast_large", "rabbit_density_mean.tif"), map_type,
                             habitat_static = file.path(ed, "lynx_habitat_corine_500_large.tif"),
                             ndvi_mean = file.path(ed, "NDVI_large", "harmonised", "mean.tif"),
                             ndvi_sd = file.path(ed, "NDVI_large", "harmonised", "sd.tif"),
                             height = file.path(ed, "vegetation_height_500_large.tif"),
                             barrier = file.path(ed, "lynx_barrier_500_large.tif"),
                             roads = file.path(ed, "lynx_roads_500_large.tif"),
                             park = file.path(ed, "park_mask_500_large.tif"),
                             presence_dir = file.path("data", "lynx", "presence"),
                             counts_file = file.path("data", "lynx", "census", "donana_breeding_females.csv"),
                             start_year = start_year, end_year = end_year)
  post <- read.csv(posteriors[[map_type]])
  t0 <- Sys.time()
  pp <- lynx_posterior_predictive(post, inp, exe, rates, n_draws = n_draws, seed = 1, fit_years = fit_years)
  cat(sprintf("%s: %d runs in %.1f min\n", map_type, pp$n_runs, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ev <- lynx_evaluate(pp, inp, fit_years)
  list(metrics = cbind(habitat_map = map_type, ev$metrics,
                       median_d_count = stats::median(pp$stats$d_count), median_d_presence = stats::median(pp$stats$d_presence),
                       habitat_share = mean(lynx_habitat_values(inp, lapply(post[, names(lynx_priors(map_type))], stats::weighted.mean, post$weights))[inp$land, ] == 2)),
       by_year = cbind(habitat_map = map_type, ev$by_year))
})

metrics <- do.call(rbind, lapply(res, `[[`, "metrics"))
by_year <- do.call(rbind, lapply(res, `[[`, "by_year"))
write.csv(metrics, file.path(out_dir, "comparison_metrics.csv"), row.names = FALSE)
write.csv(by_year, file.path(out_dir, "by_year.csv"), row.names = FALSE)
print(metrics, digits = 3)

p <- ggplot(by_year, aes(year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = habitat_map), alpha = 0.2) +
  geom_line(aes(y = median, colour = habitat_map), linewidth = 0.9) +
  geom_point(aes(y = obs, shape = period), na.rm = TRUE) +
  scale_shape_manual(values = c(calibration = 16, "held out" = 1, "spin-up" = 4)) +
  geom_vline(xintercept = max(fit_years) + 0.5, linetype = 2, colour = "grey50") +
  labs(x = NULL, y = "Territorial / breeding females", colour = "Habitat map", fill = "Habitat map", shape = "Observed",
       title = "Doñana lynx: CORINE vs NDVI habitat map",
       subtitle = sprintf("Posterior predictive median and 90%% interval (%d runs each); dashed: end of calibration", n_draws)) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "comparison.png"), p, width = 9, height = 4.5)
