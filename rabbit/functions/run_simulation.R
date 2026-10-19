

run_simulation <- function(posterior_df, start_year, end_year,
                           Cov_list, DT_template,
                           obs_transects, obs_df,
                           demograph_mod, rabbit_demographic_rates,
                           N_samples = 3,
                           compare_obs = F,
                           demographic_climate){

  required_cols <- c("K_max", "ndvi_opt", "ndvi_sigma", "ndvi_sd_opt", "ndvi_sd_sigma",
                     "height_opt", "height_sigma", "dens_opt", "R_lambda", "R_sigma",
                     "ndvi_w", "ndvi_sd_w", "height_w", "weights")
  missing_cols <- setdiff(required_cols, names(posterior_df))
  if (length(missing_cols) > 0) {
    stop("posterior_df is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Posteriors from before the suitability weights were added must be re-calibrated.")
  }

  tmpdir <- tempfile(tmpdir = tempdir())   # unique path, not yet created
  dir.create(tmpdir, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  idx <- sample(x = nrow(posterior_df), size = N_samples, prob = posterior_df$weights, replace = T)

  all_reps <- vector("list", N_samples)

  for (i in c(1:length(idx))) {

    dens_opt       <- round(posterior_df$dens_opt[idx[i]])   # discrete
    R_lambda       <- posterior_df$R_lambda[idx[i]]
    R_sigma        <- posterior_df$R_sigma[idx[i]]

    map_param_list <- list(
      K      = posterior_df$K_max[idx[i]],
      N_opt  = posterior_df$ndvi_opt[idx[i]],
      N_sig  = posterior_df$ndvi_sigma[idx[i]],
      S_opt  = posterior_df$ndvi_sd_opt[idx[i]],
      S_sig  = posterior_df$ndvi_sd_sigma[idx[i]],
      H_opt  = posterior_df$height_opt[idx[i]],
      H_sig  = posterior_df$height_sigma[idx[i]],
      N_w    = posterior_df$ndvi_w[idx[i]],
      S_w    = posterior_df$ndvi_sd_w[idx[i]],
      H_w    = posterior_df$height_w[idx[i]]
    )

    hab_folder <- file.path(tmpdir, idx[i], "habitat_folder")

    k_file <- file.path(hab_folder, "K_maps.kbin")

    if(!dir.exists(hab_folder)) {
      dir.create(hab_folder, recursive = T)

      # One binary file with the K maps of all months. Months after the last one reuse the
      # last month in the Pascal program, so no copy of December into the next January is needed.
      create_K_map_gaussian(
        Cov_list   = Cov_list,
        param_list = map_param_list,
        output_dir = hab_folder,
        start_year = start_year,
        end_year   = end_year
      )
    }

    sim_dir <- file.path(tmpdir, paste("sim", i, sep = "_"))
    dir.create(sim_dir, showWarnings = FALSE)

    cmd_args <- c(
      demograph_mod,
      sim_dir,
      as.character(start_year),
      as.character(end_year),
      rabbit_demographic_rates,
      k_file,                     # first habitat map (the Pascal program reads K from the file below)
      k_file,                     # K maps of all months
      demographic_climate,
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

    if(compare_obs){
      # get_sim_data() is sourced from function_rabbit_distances.R
      sim_data <- get_sim_data(sim_output = sim_dir,
                               end_year = end_year,
                               start_year = start_year,
                               # DetP_dir = detect_folder,
                               DT_template = DT_template,
                               transectsD = obs_transects,
                               obs_data = obs_df)
      all_reps[[i]] <- sim_data
      }


  }

  maps <- calculate_matrix_statistics(folder_paths = tmpdir,
                                      template_raster = DT_template)

  df <- if (compare_obs) {
    bind_rows(all_reps, .id = "id") %>% mutate(id = as.integer(id))
  } else {
    NULL
  }

    return(list("pop_maps" = maps,
                "sim_df" = df))
  }


