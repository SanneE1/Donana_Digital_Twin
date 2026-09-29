# ----------------------------------------------------------------------------------------
# ABC model
# ----------------------------------------------------------------------------------------

model_fn <- function(params, start_yr, end_yr) {
  
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
    H_sig     = params[8]
  )

  dens_opt   <- round(params[9])   # discrete
  R_lambda   <- params[10]
  R_sigma    <- params[11]
  
  # obs_p      <- params[8]
  # obs_p_ndvi <- params[9]
  
  tryCatch({
    tmpdir <- tempfile(tmpdir = tempdir())   # unique path, not yet created
    dir.create(tmpdir, recursive = TRUE)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    
    hab_folder <- file.path(tmpdir, "habitat_folder")
    # detect_folder <- file.path(tmpdir, "detection_folder")
    
    dir.create(hab_folder)
    # dir.create(detect_folder)
    
    create_K_map_gaussian(
      Cov_list   = Cov_list, 
      param_list = map_param_list,
      output_dir = hab_folder,
      start_year = start_yr, 
      end_year   = end_yr
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
    
    cmd_args <- c(
      demograph_mod,
      sim_dir,
      as.character(start_yr),
      as.character(end_yr),
      rabbit_demographic_rates,
      file.path(hab_folder, list.files(hab_folder)[1]),
      hab_folder,
      out_dir_clim,
      "0",                        # Do not create all monthly maps (only census months)
      "0",                        # Do not use a starting map, but use carrying capacity map to initiate
      as.character(dens_opt),
      as.character(R_lambda),
      as.character(R_sigma)
    )
    
    result <- system2(cmd_args[1], args = cmd_args[-1],
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
    
    merged <- left_join(obs_data_parsed, sim_data,
                        by = c("Fecha", "transect"),
                        suffix = c("_obs", "_sim")) %>%
      arrange(Fecha, transect)
    
    if (nrow(merged) == 0) {
      message("No matched rows for params: ", paste(round(params, 4), 
                                                    collapse = ", "))
      return(c(-5, -5))
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
      return(c(r_term, change_term))
    }else{
      return(c(-5, -5))  
    }
    
  }, error = function(e) {
    message("Failed with params: ", paste(round(params, 4), 
                                          collapse = ", "))
    message("Model run failed: ", conditionMessage(e))
    return(c(-5, -5))
  })
}
