create_K_map_gaussian <- function(Cov_list, param_list,
                                  output_dir,
                                  start_year, end_year,
                                  plot_K = F) {

  Cov_list <- lapply(Cov_list, function(x){
    if(!is.null(time(x)) && !all(is.na(time(x)))){
      time(x) <- floor_date(time(x), "month")
    }
    return(x)
  })

  
  for (year in start_year:end_year) {
    for (month in 1:12){
        idn <- format(time(Cov_list$ndvi_mean), "%Y-%m") ==
          sprintf("%04d-%02d", year, month)
        idf <- format(time(Cov_list$flood), "%Y-%m") ==
          sprintf("%04d-%02d", year, month)
        idh <- if(all(is.na(time(Cov_list$height)))){
          seq_len(nlyr(Cov_list$height)) == 1
        } else {
          format(time(Cov_list$height), "%Y-%m") ==
            sprintf("%04d-%02d", year, month)
        }

        ndvi_mean <- subset(Cov_list$ndvi_mean, idn)
        ndvi_sd   <- subset(Cov_list$ndvi_sd, idn)

        height <- subset(Cov_list$height, idh)
        flood  <- tryCatch(subset(Cov_list$flood, idf),
                           error = function(e) NULL)

        # Each covariate contributes a suitability index in [0, 1], peaking at
        # its optimum; K is their weighted average scaled by K_max, so the
        # result is always bounded in [0, K_max] regardless of the sampled
        # sigmas. Weights are normalised to sum to 1.
        suit_ndvi   <- exp(-(ndvi_mean - param_list$N_opt)^2 / (2 * param_list$N_sig^2))
        suit_sd     <- exp(-(ndvi_sd   - param_list$S_opt)^2 / (2 * param_list$S_sig^2))
        suit_height <- exp(-(height    - param_list$H_opt)^2 / (2 * param_list$H_sig^2))

        w <- c(param_list$N_w, param_list$S_w, param_list$H_w)
        w <- w / sum(w)

        Krast <- param_list$K * (w[1] * suit_ndvi + w[2] * suit_sd + w[3] * suit_height) + 1

        if(!is.null(flood)){
          Krast <- ifel(is.na(flood), Krast, ifel(flood == 1, 0, Krast))
        }

        Krast[is.na(Krast)] <- 0 
        
        if(plot_K) {
          plot(Krast, main = sprintf("%04d-%02d", year, month))
        }
        
        save_raster_as_input_format(Krast, file.path(output_dir,
                                                     paste0(year, "_", sprintf("%02d", month), ".txt")))

      }
  }
}
