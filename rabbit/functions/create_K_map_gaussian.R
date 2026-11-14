# Carrying capacity K per cell and month from Gaussian suitability of NDVI mean, NDVI SD and
# vegetation height. Needs K_map_binary.R (prepare_K_covariates, write_K_binary, write_K_text).
#
# format = "binary": writes one file, <output_dir>/K_maps.kbin, and returns its path. This is
#                    what the Pascal program reads when its habitat argument is a file.
# format = "text":   writes one <year>_<month>.txt per month (older input, easier to inspect)
#                    and returns output_dir.
# Cov_list can also be the result of prepare_K_covariates(), which saves repeating that step
# when many parameter sets are run on the same covariates.
# study_area: raster on the same grid as the covariates (e.g. the template), NA outside the study
# area. When given, any cell inside it without NDVI mean/SD in a month is an error: such cells
# would otherwise silently get K = 0 (fill the gaps first, see build_monthly_ndvi()).
create_K_map_gaussian <- function(Cov_list, param_list,
                                  output_dir,
                                  start_year, end_year,
                                  study_area = NULL,
                                  plot_K = F,
                                  format = c("binary", "text")) {

  format <- match.arg(format)

  cov <- if (is.null(Cov_list$months)) prepare_K_covariates(Cov_list, start_year, end_year) else Cov_list

  months_missing <- which(colSums(!is.na(cov$ndvi_mean)) == 0 | colSums(!is.na(cov$ndvi_sd)) == 0)
  if (length(months_missing) > 0) {
    stop("No NDVI mean/SD layer for month(s): ", paste(cov$months[months_missing], collapse = ", "))
  }

  if (!is.null(study_area)) {
    inside <- !is.na(terra::values(study_area[[1]], mat = FALSE))
    if (length(inside) != nrow(cov$ndvi_mean)) {
      stop("study_area does not have the same number of cells as the covariates")
    }
    n_gaps <- colSums(is.na(cov$ndvi_mean[inside, , drop = FALSE]) | is.na(cov$ndvi_sd[inside, , drop = FALSE]))
    if (any(n_gaps > 0)) {
      gaps <- which(n_gaps > 0)
      stop("Cells inside the study area without NDVI mean/SD (they would get K = 0) in ",
           length(gaps), " month(s): ",
           paste0(cov$months[gaps], " (", n_gaps[gaps], " cells)", collapse = ", "),
           ". Fill the gaps first, e.g. with build_monthly_ndvi().")
    }
  }

  # Each covariate contributes a suitability index in [0, 1], peaking at its optimum;
  # K is their weighted average scaled by K_max, so the result is always bounded in
  # [0, K_max] regardless of the sampled sigmas. Weights are normalised to sum to 1.
  w <- c(param_list$N_w, param_list$S_w, param_list$H_w)
  w <- w / sum(w)

  K <- param_list$K * (
    w[1] * exp(-(cov$ndvi_mean - param_list$N_opt)^2 / (2 * param_list$N_sig^2)) +
    w[2] * exp(-(cov$ndvi_sd   - param_list$S_opt)^2 / (2 * param_list$S_sig^2)) +
    w[3] * exp(-(cov$height    - param_list$H_opt)^2 / (2 * param_list$H_sig^2))
  ) + 1

  # Flooded cells are not habitat; where the flood map has no value the cell is kept
  K[!is.na(cov$flood) & cov$flood == 1] <- 0
  K[is.na(K)] <- 0
  # Truncate, as the earlier writeRaster(datatype = "INT2S") route did, so K is the same as in
  # calibrations run before the binary format.
  K <- floor(K)

  if (plot_K) {
    for (i in seq_len(ncol(K))) {
      r <- cov$template
      terra::values(r) <- K[, i]
      plot(r, main = cov$months[i])
    }
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (format == "binary") {
    return(write_K_binary(K, ncol = cov$ncol, nrow = cov$nrow,
                          first_year = start_year, first_month = 1,
                          file = file.path(output_dir, "K_maps.kbin")))
  }

  for (i in seq_len(ncol(K))) {
    write_K_text(K[, i], cov$ncol, cov$nrow,
                 file.path(output_dir, paste0(sub("-", "_", cov$months[i]), ".txt")))
  }
  invisible(output_dir)
}
