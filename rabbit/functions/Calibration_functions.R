# ----------------------------------------------------------------------------------------
# ABC model
# ----------------------------------------------------------------------------------------

# The model is called with named arguments (rabbit_model_args.R).
# params[1] is the seed EasyABC passes; it is also the model's random seed, so a
# simulation can be repeated.
# RHD (optional, see rhd_settings.R and rhd_summaries.R):
#   rhd_settings     named list of RHD settings for the Pascal program, including RHD_mode
#                    (NULL = those of the rates file)
#   rhd_param_names  names of RHD settings calibrated here; their values are params[15], params[16], ...
#                    "log10_<name>" is calibrated on the log10 scale (e.g. log10_RHD_spark);
#                    "RHDMort" sets RHDMort_L = RHDMort_U (schedule mode). See model_param_settings().
#   rhd_target       list(scale = "park" or "cell", target = c(24, 36), from = month index or NULL);
#                    when given, a third statistic is returned: the distance of the outbreak
#                    intervals from the target range (rhd_interval_distance), 0 = inside it
#   return_sim       TRUE: return the simulated and observed KAI per census (for model comparison)
#                    instead of the summary statistics
model_fn <- function(params, start_yr, end_yr,
                     Cov_list, DT_template,
                     demograph_mod, rabbit_demographic_rates,
                     obs_transects, obs_df, out_dir_clim,
                     rhd_settings = NULL, rhd_param_names = NULL, rhd_target = NULL,
                     return_sim = FALSE) {

  # params = c(1, 1.1,1.1,20,1,0.01,0.1,0,0)   # just here for quick test runs
  set.seed(params[1])
  # alpha      <- params[2]
  # beta_param <- params[3]
  map_param_list <- list(
    K         = round(params[2]),   # discrete; prior samples a continuous uniform
    N_opt     = params[3],
    N_sig     = params[4],
    S_opt     = params[5],
    S_sig     = params[6],
    H_opt     = params[7],
    H_sig     = params[8],
    N_w       = params[12],         # suitability weights (normalised in create_K_map_gaussian)
    S_w       = params[13],
    H_w       = params[14]
  )

  dens_opt   <- round(params[9])   # discrete
  R_lambda   <- params[10]
  R_sigma    <- params[11]
  
  # obs_p      <- params[8]
  # obs_p_ndvi <- params[9]
  
  if (length(rhd_param_names) > 0) {
    rhd_values <- stats::setNames(as.list(params[14 + seq_along(rhd_param_names)]), rhd_param_names)
    rhd_settings <- utils::modifyList(as.list(rhd_settings), model_param_settings(rhd_values))
  }
  failed <- if (is.null(rhd_target)) c(-5, -5) else c(-5, -5, 5)
  
  tryCatch({
    tmpdir <- tempfile(tmpdir = tempdir())   # unique path, not yet created
    dir.create(tmpdir, recursive = TRUE)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    
    hab_folder <- file.path(tmpdir, "habitat_folder")
    # detect_folder <- file.path(tmpdir, "detection_folder")
    
    dir.create(hab_folder)
    # dir.create(detect_folder)
    
    # One binary file with the K maps of all months (see K_map_binary.R)
    k_file <- create_K_map_gaussian(
      Cov_list   = Cov_list,
      param_list = map_param_list,
      output_dir = hab_folder,
      start_year = start_yr,
      end_year   = end_yr,
      study_area = DT_template
      )
    
    # create_detectionP_map(
    #   Cov_list   = Cov_list,
    #   output_dir = detect_folder,
    #   obs_p = obs_p,
    #   obs_p_ndvi = obs_p_ndvi,
    #   start_year = start_yr,
    #   end_year   = end_yr
    # )
    
    sim_dir <- file.path(tmpdir, "sim")
    dir.create(sim_dir, showWarnings = FALSE)
    
    # Named arguments: run settings, then the dispersal and RHD settings of this
    # parameter set (they override the rates file). Only census-month maps; the
    # population starts at carrying capacity.
    args <- rabbit_model_args(c(
      list(output_dir     = sim_dir,
           start_year     = start_yr,
           end_year       = end_yr,
           rates_file     = rabbit_demographic_rates,
           habitat_folder = k_file,
           climate_folder = out_dir_clim,
           seed           = as.integer(params[1]),
           dens_opt       = dens_opt,
           lambda         = R_lambda,
           sigma          = R_sigma),
      rhd_settings))

    result <- system2(demograph_mod, args = args,
                      stdout = TRUE, stderr = TRUE,
                      timeout = 600)
    
    if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
      stop("Subprocess failed: ", paste(result, collapse = "\n"))
    }
    
    # get_sim_data() is sourced from function_rabbit_distances.R
    sim_data <- get_sim_data(sim_output = sim_dir, 
                             start_year = start_yr, 
                             end_year = end_yr, 
                             # DetP_dir = detect_folder,
                             DT_template = DT_template, 
                             transectsD = obs_transects, 
                             obs_data = obs_df)
    
    
    obs_data_parsed <- obs_df %>%
      dplyr::select(Fecha, KAI, transect) %>%
      mutate(Fecha = dmy(Fecha))
    
    merged <- left_join(sim_data, obs_data_parsed,
                        by = c("Fecha", "transect"),
                        suffix = c("_sim", "_obs")) %>%
      arrange(Fecha, transect)

    if (return_sim) return(merged)

    if (nrow(merged) == 0 || all(is.na(merged$KAI_obs))) {
      message("No matched rows for params: ", paste(round(params, 4), 
                                                    collapse = ", "))
      return(failed)
    }
    
    suppressWarnings({
      metrics <- merged %>%
        group_by(transect) %>%
        arrange(Fecha, .by_group = TRUE) %>%   # ensure chronological order
        summarise(
          rmse = sqrt(mean((KAI_obs - KAI_sim)^2, na.rm = T)),    # overall error
          rmse_scaled = rmse / sd(KAI_obs, na.rm = T),            # scaled data fit
          r = cor(KAI_obs, KAI_sim, use = "complete.obs"),        # overal fit/pattern match
          r_delta = cor(diff(KAI_obs), diff(KAI_sim),             # Year-to-year change fidelity
                        use = "complete.obs")         
        )
    })
    
    
    rmse_term <- mean(metrics$rmse_scaled, na.rm = TRUE)
    r_term <- median(metrics$r, na.rm = TRUE)
    change_term <- median(metrics$r_delta, na.rm = TRUE)
    
    if(is.finite(rmse_term) & is.finite(r_term) & is.finite(change_term)){
      if (is.null(rhd_target)) return(c(r_term, change_term))
      rhd_intervals_sim <- rhd_intervals(read_rhd_logs(sim_dir), scale = rhd_target$scale,
                                         from = rhd_target$from)
      return(c(r_term, change_term, rhd_interval_distance(rhd_intervals_sim, rhd_target$target)))
    }else{
      return(failed)
    }
    
  }, error = function(e) {
    message("Failed with params: ", paste(round(params, 4), 
                                          collapse = ", "))
    message("Model run failed: ", conditionMessage(e))
    return(failed)
  })
}
