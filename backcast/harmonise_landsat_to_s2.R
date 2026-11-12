#
# Backcast phase 2: adjust Landsat NDVI to Sentinel-2
#
# The rabbit model was calibrated on the monthly 500 m mean and SD of 10 m Sentinel-2 NDVI. Landsat
# (30 m, other bands and sensors) gives different values for the same surface, so each Landsat
# sensor is mapped onto Sentinel-2 with a linear model per statistic, fitted on matching cells and
# months where both are available:
#   Landsat 8 -> Sentinel-2   overlap 2015-2025
#   Landsat 7 -> Sentinel-2   overlap 2015-2022
#   Landsat 5 -> adjusted Landsat 7   overlap 1999-2011 (Landsat 5 stopped before Sentinel-2)
# Each fit is first checked on held-out years, then refitted on all overlapping months. Fits use
# observed values only: the raw monthly products, without cells with too few clear pixels and
# without gap filling.
#
# Output (data/backcast):
#   monthly_ndvi_harmonised/mean.tif, sd.tif   adjusted Landsat before Sentinel-2, Sentinel-2 after,
#                                              gap-filled as in build_monthly_ndvi()
#   monthly_ndvi_harmonised/filled.tif         1 where filled, 0 where observed
#   harmonisation_fits.csv                     coefficients and hold-out error per sensor/statistic
#   harmonisation_timeseries.png               study-area mean of each series over time
#
# Run environmental_data/rebuild_monthly_ndvi.R and Workflow_backcast_data.R (phase = 2) first.

library(terra)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

bc_dir      <- file.path("data", "backcast")

source(file.path("backcast", "functions", "backcast_helpers.R"))
source(file.path("environmental_data", "functions", "build_monthly_ndvi.R"))

# Area: "large" (the default: lynx area, template_raster_500_large.tif, products of
# build_monthly_ndvi_500m() from the Google Earth Engine downloads) or "park" (rabbit grid). Set
# harmonise_area before sourcing this script to choose.
# The park Sentinel-2 product (monthly_ndvi_v2) was removed on 2026-11-10: the park now uses
# NDVI_large/harmonised, cut to the park template, so "park" only works if that product is rebuilt.
if(!exists("harmonise_area")) harmonise_area <- "large"
if(harmonise_area == "park"){
  if(!dir.exists(file.path("data", "environmental_data", "monthly_ndvi_v2")))
    stop("harmonise_area = 'park' needs data/environmental_data/monthly_ndvi_v2, which was removed; use 'large'")
  s2_dir      <- file.path("data", "environmental_data", "monthly_ndvi_v2")
  landsat_dir <- function(p) file.path(bc_dir, paste0("monthly_ndvi_", p))
  out_dir     <- file.path(bc_dir, "monthly_ndvi_harmonised")
  fits_file   <- file.path(bc_dir, "harmonisation_fits.csv")
  plot_file   <- file.path(bc_dir, "harmonisation_timeseries.png")
} else if(harmonise_area == "large"){
  large_dir   <- file.path("data", "environmental_data", "NDVI_large")
  s2_dir      <- file.path(large_dir, "sentinel-2")
  landsat_dir <- function(p) file.path(large_dir, p)
  # harmonise_out_dir (set before sourcing) writes elsewhere, e.g. harmonised_1984 for the backcast
  out_dir     <- if(exists("harmonise_out_dir")) harmonise_out_dir else file.path(large_dir, "harmonised")
  fits_file   <- file.path(out_dir, "harmonisation_fits.csv")
  plot_file   <- file.path(out_dir, "harmonisation_timeseries.png")
} else stop("harmonise_area must be \"park\" or \"large\"")
min_valid_frac <- 0.25
dir.create(out_dir, showWarnings = FALSE)

# Observed values of a monthly product: raw stack without cells below min_valid_frac
read_observed <- function(product_dir, stat){
  r <- read_monthly_ndvi(file.path(product_dir, paste0("raw_", stat, ".tif")))
  valid <- rast(file.path(product_dir, "valid_fraction.tif"))
  r[valid < min_valid_frac] <- NA
  r
}

first_s2 <- as.Date("2015-07-01")

# Hold-out years per sensor, for the error reported in harmonisation_fits.csv
holdout <- list("landsat-8" = 2023:2025, "landsat-7" = 2021:2022, "landsat-5" = 2010:2011)

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------

# cell x month table of two monthly stacks, for the months both have
paired_values <- function(x, y){
  tx <- as.Date(time(x)); ty <- as.Date(time(y))
  common <- sort(as.Date(intersect(tx, ty)))
  if(length(common) == 0) stop("No overlapping months")
  vx <- values(x[[match(common, tx)]], mat = TRUE)
  vy <- values(y[[match(common, ty)]], mat = TRUE)
  data.frame(cell  = rep(seq_len(nrow(vx)), times = length(common)),
             month = rep(common, each = nrow(vx)),
             x = as.vector(vx), y = as.vector(vy)) %>%
    filter(!is.na(x), !is.na(y))
}

fit_sensor <- function(x, y, sensor, stat, holdout_years){
  d <- paired_values(x, y)
  train <- d %>% filter(!(year(month) %in% holdout_years))
  test  <- d %>% filter(year(month) %in% holdout_years)

  m_train <- lm(y ~ x, data = train)
  pred    <- predict(m_train, newdata = test)
  m_all   <- lm(y ~ x, data = d)

  list(model = m_all,
       summary = data.frame(
         sensor = sensor, stat = stat,
         intercept = coef(m_all)[1], slope = coef(m_all)[2],
         r2 = summary(m_all)$r.squared,
         n_pairs = nrow(d), n_months = length(unique(d$month)),
         holdout_years = paste(range(holdout_years), collapse = "-"),
         holdout_n = nrow(test),
         holdout_rmse = sqrt(mean((pred - test$y)^2)),
         holdout_bias = mean(pred - test$y),
         holdout_rmse_unadjusted = sqrt(mean((test$x - test$y)^2)),
         row.names = NULL))
}

apply_fit <- function(r, model, stat){
  out <- coef(model)[1] + coef(model)[2] * r
  out <- if(stat == "mean") clamp(out, -1, 1) else clamp(out, 0, Inf)
  time(out) <- time(r)
  out
}

#-------------------------------------------------------------------------------
# Fit and apply per statistic
#-------------------------------------------------------------------------------

fits <- list()
series <- list()

for(stat in c("mean", "sd")){
  s2  <- read_observed(s2_dir, stat)
  ls  <- lapply(c(l5 = "landsat-5", l7 = "landsat-7", l8 = "landsat-8"), function(p){
    read_observed(landsat_dir(p), stat)
  })

  f8 <- fit_sensor(ls$l8, s2, "landsat-8", stat, holdout$`landsat-8`)
  f7 <- fit_sensor(ls$l7, s2, "landsat-7", stat, holdout$`landsat-7`)
  l7_adj <- apply_fit(ls$l7, f7$model, stat)
  f5 <- fit_sensor(ls$l5, l7_adj, "landsat-5 (via landsat-7)", stat, holdout$`landsat-5`)

  adj <- list(l8 = apply_fit(ls$l8, f8$model, stat),
              l5 = apply_fit(ls$l5, f5$model, stat),
              l7 = l7_adj)
  fits[[stat]] <- bind_rows(f8$summary, f7$summary, f5$summary)

  # Months before Sentinel-2: per cell the first available of Landsat 8, 5, 7. Landsat 7 comes
  # last because of the striped gaps in its images after May 2003.
  pre_months <- sort(unique(do.call(c, lapply(adj, function(r) as.Date(time(r))))))
  pre_months <- pre_months[pre_months < first_s2]

  pre <- rast(lapply(pre_months, function(m){
    layers <- lapply(adj, function(r){
      i <- which(as.Date(time(r)) == m)
      if(length(i) == 1) r[[i]] else NULL
    })
    Reduce(cover, Filter(Negate(is.null), layers))
  }))
  time(pre) <- pre_months

  s2_keep <- s2[[as.Date(time(s2)) >= first_s2]]
  out <- c(pre, s2_keep)
  time(out) <- c(pre_months, as.Date(time(s2_keep)))
  
  g <- fill_ndvi_gaps(out)
  writeRaster(g$filled, file.path(out_dir, paste0(stat, ".tif")), overwrite = TRUE)
  if(stat == "mean") writeRaster(g$missing * 1, file.path(out_dir, "filled.tif"), overwrite = TRUE, datatype = "INT1U")

  # study-area means for the plot
  area_mean <- function(r, name){
    data.frame(month = as.Date(time(r)), value = global(r, "mean", na.rm = TRUE)[, 1], series = name)
  }
  series[[stat]] <- bind_rows(area_mean(s2, "Sentinel-2"),
                              area_mean(ls$l8, "Landsat 8 (raw)"), area_mean(adj$l8, "Landsat 8 (adjusted)"),
                              area_mean(ls$l7, "Landsat 7 (raw)"), area_mean(adj$l7, "Landsat 7 (adjusted)"),
                              area_mean(ls$l5, "Landsat 5 (raw)"), area_mean(adj$l5, "Landsat 5 (adjusted)")) %>%
    mutate(stat = paste("NDVI", stat))
}

fit_table <- bind_rows(fits)
write.csv(fit_table, fits_file, row.names = FALSE)
print(fit_table)

plot_df <- bind_rows(series) %>%
  mutate(sensor = sub(" \\(.*", "", series),
         type = ifelse(grepl("raw", series), "raw", "adjusted / Sentinel-2"))

p <- ggplot(plot_df, aes(month, value, colour = sensor, linetype = type)) +
  geom_line(linewidth = 0.3) +
  geom_vline(xintercept = first_s2, linetype = "dotted", colour = "grey50") +
  facet_wrap(vars(stat), ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Study-area mean", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11)

ggsave(p, filename = plot_file, width = 10, height = 6)
