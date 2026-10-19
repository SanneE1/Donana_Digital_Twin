# ----------------------------------------------------------------------------------------
# Carrying capacity (K) maps for the Pascal rabbit model, stored in ONE binary file.
#
# Any way of calculating K can use this: build a cells x months matrix of integer K values
# (cells in raster order: row 1 first, left to right; months consecutive from first_year /
# first_month) and call write_K_binary(). The Pascal program only reads this file, so it does
# not care how K was calculated.
#
# File layout (little endian):
#   bytes  0-3   "KMAP"
#   int32 x 7    version (1), ncol, nrow, n_months, first_year, first_month, bytes_per_value (1 or 2)
#   then n_months blocks, each nrow * ncol values (row 1 = top row, x fastest), unsigned 8 bit
#   when K <= 255 everywhere, otherwise unsigned 16 bit
# Months after the last one stored reuse the last month in the Pascal program.
# ----------------------------------------------------------------------------------------

write_K_binary <- function(K_mat, ncol, nrow, first_year, first_month = 1, file) {
  K_mat <- as.matrix(K_mat)
  if (nrow(K_mat) != ncol * nrow) {
    stop("K_mat has ", nrow(K_mat), " rows but ncol * nrow = ", ncol * nrow)
  }
  if (anyNA(K_mat)) stop("K_mat contains NA; set non-habitat cells to 0 first")
  if (min(K_mat) < 0) stop("K_mat contains negative values; non-habitat cells must be 0")

  K_mat <- round(K_mat)
  bytes <- if (max(K_mat) <= 255) 1L else 2L
  if (max(K_mat) > 65535) stop("K values above 65535 are not supported")

  con <- file(file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("KMAP"), con)
  writeBin(as.integer(c(1L, ncol, nrow, ncol(K_mat), first_year, first_month, bytes)),
           con, size = 4L, endian = "little")
  if (bytes == 1L) {
    writeBin(as.raw(as.vector(K_mat)), con)
  } else {
    writeBin(as.integer(as.vector(K_mat)), con, size = 2L, endian = "little")
  }
  invisible(file)
}


# One month as the text format the Pascal program also accepts (ncol nrow header, then rows).
write_K_text <- function(K_vec, ncol, nrow, file) {
  m <- matrix(as.integer(round(K_vec)), nrow = nrow, ncol = ncol, byrow = TRUE)
  lines <- c(paste(ncol, nrow), do.call(paste, c(as.data.frame(m), sep = " ")))
  writeLines(lines, file)
  invisible(file)
}


# Covariate layers as plain cells x months matrices for each target month, so K can be
# calculated with array maths. Covariates that are static (one layer, or without dates) are
# repeated for every month. Missing layers become NA columns.
prepare_K_covariates <- function(Cov_list, start_year, end_year) {
  months <- format(seq(as.Date(sprintf("%04d-01-01", start_year)),
                       as.Date(sprintf("%04d-12-01", end_year)), by = "month"), "%Y-%m")

  as_month_matrix <- function(r, allow_static) {
    nc <- terra::ncell(r)
    tm <- terra::time(r)
    if (allow_static && (is.null(tm) || all(is.na(tm)))) {
      v <- terra::values(r[[1]], mat = FALSE)
      return(matrix(v, nrow = nc, ncol = length(months)))
    }
    idx <- match(months, format(lubridate::floor_date(tm, "month"), "%Y-%m"))
    out <- matrix(NA_real_, nrow = nc, ncol = length(months))
    ok <- !is.na(idx)
    if (any(ok)) out[, ok] <- terra::values(r[[idx[ok]]], mat = TRUE)
    out
  }

  list(ndvi_mean = as_month_matrix(Cov_list$ndvi_mean, FALSE),
       ndvi_sd   = as_month_matrix(Cov_list$ndvi_sd, FALSE),
       height    = as_month_matrix(Cov_list$height, TRUE),
       flood     = as_month_matrix(Cov_list$flood, FALSE),
       months    = months,
       ncol      = terra::ncol(Cov_list$ndvi_mean),
       nrow      = terra::nrow(Cov_list$ndvi_mean),
       template  = Cov_list$ndvi_mean[[1]])
}
