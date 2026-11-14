# ----------------------------------------------------------------------------
# Running the lynx IBM (lynx/Pascal_program, Model_cmd) from R
#
# The program is called as   <exe> <settings file> <output dir>
# with a settings file of "key value" lines (start_year, end_year, seed, lynx_demography,
# map_lynx_pops, lynx_start_size, lynx_reintro, habitat_folder, breeding_file / breeding_folder,
# prey_folder, all_year_maps, create_additional_output, ...). It writes lynx_territories.csv
# (lynx per population ID, sex and status in May of every year), lynx_pop_size.csv and, with
# all_year_maps 1, maps/FemalesMap_status_yr_<year>.csv and MalesMap_status_yr_<year>.csv
# (-1 = no territory, 2 = settling, 3 = territorial; rows north to south).
# ----------------------------------------------------------------------------

# Writes cell values (template order: rows north to south) in the IBM map format
write_lynx_values <- function(v, ncol, nrow, file){
  v[is.na(v)] <- -9999
  m <- matrix(as.integer(round(v)), nrow = nrow, byrow = TRUE)
  lines <- c(paste(ncol, nrow), apply(m, 1, paste, collapse = " "))
  writeLines(lines, file)
  invisible(file)
}

# Writes a settings file; `settings` is a named list
write_lynx_settings <- function(settings, file){
  vals <- vapply(settings, function(v) if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE) else as.character(v), "")
  writeLines(paste(names(settings), vals), file)
  invisible(file)
}

# Demographic rates file: the base file (key value lines, surv_<age> with 6 columns) with the
# calibrated values put in.
#   rates       named list of single values replacing keys of the base file (e.g. Tsize, rep_prob,
#               surv_disp_rho)
#   surv_adult  named numeric vector: annual survival from age 2 on, per survival column
#               (c("0" = outside, "2" = Doñana)); ages above max_age keep 0
read_lynx_rates <- function(file){
  l <- trimws(sub("^﻿", "", readLines(file, encoding = "UTF-8", warn = FALSE)))
  l <- l[nzchar(l)]
  sp <- strsplit(l, "\\s+")
  stats::setNames(lapply(sp, function(x) x[-1]), vapply(sp, `[`, "", 1))
}

write_lynx_rates <- function(base, file, rates = list(), surv_adult = NULL){
  if (is.character(base)) base <- read_lynx_rates(base)
  for (k in names(rates)) {
    if (is.null(base[[k]])) stop("Unknown lynx rate: ", k)
    base[[k]] <- format(rates[[k]], scientific = FALSE, trim = TRUE)
  }
  if (!is.null(surv_adult)) {
    max_age <- as.integer(base$max_age)
    for (a in 2:21) {
      key <- paste0("surv_", a)
      s <- as.numeric(base[[key]])
      for (col in names(surv_adult)) s[as.integer(col) + 1] <- if (a <= max_age) surv_adult[[col]] else 0
      base[[key]] <- format(s, scientific = FALSE, trim = TRUE)
    }
  }
  writeLines(paste(names(base), vapply(base, paste, "", collapse = " ")), file)
  invisible(file)
}

# Start population file: n lynx in n_sites groups. `cells` (cell numbers of the template, e.g.
# habitat inside the observed presence area) are split into n_sites spatial clusters (k-means on
# the coordinates, fixed start) and each group starts in the cell of its cluster nearest the
# cluster centre. The IBM gives them age 3-5 and a random sex and lets them settle in n_cycles
# dispersal cycles before the first year. Lynx spread one per cell rarely settle: a male needs
# about 7 cells of female territories (Tsize x male_T_multiplier, 20% gaps allowed).
lynx_start_sites <- function(cells, template, n_sites = 3){
  xy <- terra::xyFromCell(template, cells)
  if (length(cells) <= n_sites) return(cells)
  set.seed(1)
  km <- stats::kmeans(xy, n_sites, nstart = 5)
  vapply(seq_len(n_sites), function(k){
    i <- which(km$cluster == k)
    cells[i[which.min(colSums((t(xy[i, , drop = FALSE]) - km$centers[k, ])^2))]]
  }, 0)
}

write_lynx_start_pop <- function(n, sites, template, file){
  n_site <- diff(round(seq(0, n, length.out = length(sites) + 1)))
  rc <- terra::rowColFromCell(template, sites)
  writeLines(c("N X Y pop", paste(n_site, rc[, 2], rc[, 1], "start")[n_site > 0]), file)
  invisible(file)
}

# Runs the IBM; returns the output directory. Stops with the program's output if it fails.
run_lynx_model <- function(exe, settings_file, out_dir, timeout = 1800){
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  res <- suppressWarnings(system2(exe, c(shQuote(settings_file), shQuote(out_dir)),
                                  stdout = TRUE, stderr = TRUE, timeout = timeout))
  status <- attr(res, "status")
  if ((!is.null(status) && status != 0) || !file.exists(file.path(out_dir, "lynx_territories.csv"))) {
    stop("Lynx model failed:\n", paste(utils::tail(res, 20), collapse = "\n"))
  }
  out_dir
}

# Territory counts per year (all population IDs summed, or per ID with by_pop = TRUE). The IBM
# only writes population IDs with lynx, so years after an extinction have no rows: with `years`
# given, every year is returned, with 0 where there were no lynx.
read_lynx_territories <- function(out_dir, by_pop = FALSE, years = NULL){
  d <- utils::read.csv(file.path(out_dir, "lynx_territories.csv"))
  if (by_pop) return(d)
  cols <- c("females_territorial", "males_territorial", "females_settling", "males_settling", "all_lynx")
  a <- if (nrow(d)) stats::aggregate(d[, cols], list(year = d$year), sum) else
    data.frame(year = integer(0), matrix(numeric(0), 0, length(cols), dimnames = list(NULL, cols)))
  if (is.null(years)) return(a)
  out <- data.frame(year = years, matrix(0, length(years), length(cols), dimnames = list(NULL, cols)))
  out[match(a$year, years, nomatch = 0), cols] <- a[a$year %in% years, cols]
  out
}

# Occupied cells per year: cell values of the template, 1 where a female or male territory
# (status >= 2) lies in May of that year, else 0 (NA outside the template). Matrix cells x years.
read_lynx_occupancy <- function(out_dir, years, template){
  land <- !is.na(terra::values(template[[1]], mat = FALSE))
  out <- vapply(years, function(y){
    occ <- rep(0, terra::ncell(template))
    for (sex in c("Females", "Males")) {
      f <- file.path(out_dir, "maps", sprintf("%sMap_status_yr_%d.csv", sex, y))
      if (!file.exists(f)) stop("No territory map for ", y, ": run with all_year_maps 1")
      m <- as.matrix(utils::read.csv(f, header = FALSE))
      occ[as.vector(t(m)) >= 2] <- 1
    }
    occ[!land] <- NA
    occ
  }, numeric(terra::ncell(template)))
  colnames(out) <- years
  out
}
