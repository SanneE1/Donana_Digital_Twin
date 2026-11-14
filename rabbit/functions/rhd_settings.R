
# RHD settings for the Pascal program (command line argument 14). See docs/model_settings.md for
# what each setting means. Settings left out keep the program's default.

# Defaults of the immunity-based disease process. Arrival year 0 = before the simulation starts,
# 9999 = never.
rhd_default_settings <- function(){
  list(
    RHD_mode            = "immunity",
    RHD_min_age         = 3,
    RHD_f_crit          = 0.5,
    RHD_f_slope         = 20,
    RHD_spark           = 0.002,
    RHD_spread          = 0.2,
    RHD_mortality       = 0.6,
    RHD2_mortality      = -1,
    RHD1_arrival_year   = 0,
    RHD1_arrival_month  = 1,
    RHD2_arrival_year   = 0,
    RHD2_arrival_month  = 1,
    RHD_seed_share      = 0.05,
    RHD2_cross_immunity = 0,
    RHD_init_immune_lo  = 0.2,
    RHD_init_immune_hi  = 0.8,
    RHD_log             = 1
  )
}

# Writes the settings as "key value" lines. `settings` is a named list, e.g. a modified
# rhd_default_settings(); numbers are written with "." as decimal mark.
write_rhd_settings <- function(settings, file){
  if (is.null(names(settings)) || any(names(settings) == "")) stop("settings must be a named list")
  values <- vapply(settings, function(v){
    if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE, digits = 15) else as.character(v)
  }, character(1))
  writeLines(paste(names(settings), values), file)
  invisible(file)
}

# Sets calibrated RHD settings: `names` are setting names, or "log10_" + a setting name for a
# setting calibrated on the log10 scale (e.g. log10_RHD_spark = -2.5 gives RHD_spark = 0.0032).
set_rhd_params <- function(settings, names, values){
  for (j in seq_along(names)) {
    if (startsWith(names[j], "log10_")) settings[[sub("^log10_", "", names[j])]] <- 10^values[j]
    else settings[[names[j]]] <- values[j]
  }
  settings
}

# The calibrated RHD settings among the columns of a posterior data frame
rhd_param_columns <- function(posterior_df){
  grep("^(log10_)?RHD", names(posterior_df), value = TRUE)
}
