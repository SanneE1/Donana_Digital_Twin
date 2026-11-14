library(terra)
library(lubridate)

# Combines the individual LAST flood maps into one map per month (maximum flood extent of all
# maps taken that month) on the template grid. The date of each map is read from the 8 digit
# date (YYYYMMDD) in its file name. Months that already have a map in out_dir are skipped unless
# overwrite = TRUE.
make_monthly_floodmaps <- function(fld_files, template, out_dir,
                                   min_year = 1984, overwrite = FALSE){

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  if(is.character(template)) template <- rast(template)

  dates <- as.Date(stringr::str_extract(basename(fld_files), "\\d{8}"), format = "%Y%m%d")
  keep  <- !is.na(dates) & year(dates) >= min_year
  fld_files <- fld_files[keep]
  dates     <- dates[keep]

  periods <- sort(unique(floor_date(dates, "month")))

  for(i in seq_along(periods)){
    p <- periods[i]
    out_file <- file.path(out_dir, paste0(year(p), "_", sprintf("%02d", month(p)), ".tif"))
    if(file.exists(out_file) && !overwrite) next

    try({
      r_files <- fld_files[floor_date(dates, "month") == p]

      r_list <- lapply(r_files, function(f) {
        r <- rast(f)
        r <- project(r, crs(template))
        r <- resample(r, template, method = "max")
        return(r)
      })

      r_max <- round(max(rast(r_list), na.rm = TRUE), digits = 0)

      writeRaster(r_max, out_file, overwrite = T, NAflag = 0)
      cat('Written', basename(out_file), 'from', length(r_files), 'map(s)', fill = T)
    })
  }

  invisible(out_dir)
}
