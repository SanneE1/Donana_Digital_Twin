
# Summaries of the RHD logs the Pascal program writes in immunity mode:
#   rhd_outbreaks.csv  one row per outbreak: year, month, x, y, strain, n_age, n_susceptible, deaths
#   rhd_monthly.csv    one row per month: year, month, strain, occupied_cells, outbreak_cells,
#                      rabbits, immune, rhd_deaths

read_rhd_logs <- function(sim_dir){
  list(outbreaks = read.csv(file.path(sim_dir, "rhd_outbreaks.csv")),
       monthly   = read.csv(file.path(sim_dir, "rhd_monthly.csv")))
}

# Months since year 0, to take differences between dates
rhd_month_index <- function(year, month) year * 12 + month

# Park-wide outbreak waves: peaks in the monthly share of occupied cells with an outbreak.
# The share is smoothed with a centred moving average of `smooth` months; a wave is a local
# maximum above `threshold` (default: the `threshold_quantile` quantile of the smoothed share);
# peaks closer than `min_gap` months are merged, keeping the highest. Returns the peak months and
# the intervals between them (months).
rhd_park_waves <- function(monthly, smooth = 3, threshold = NULL, threshold_quantile = 0.75,
                           min_gap = 6, from = NULL){
  m <- monthly[order(monthly$year, monthly$month), ]
  m$t <- rhd_month_index(m$year, m$month)
  if (!is.null(from)) m <- m[m$t >= from, ]
  share <- ifelse(m$occupied_cells > 0, m$outbreak_cells / m$occupied_cells, 0)
  s <- as.numeric(stats::filter(share, rep(1 / smooth, smooth), sides = 2))
  s[is.na(s)] <- share[is.na(s)]
  if (is.null(threshold)) threshold <- stats::quantile(s, threshold_quantile, names = FALSE)

  n <- length(s)
  is_peak <- s > threshold & s >= c(-Inf, s[-n]) & s >= c(s[-1], -Inf)
  peaks <- which(is_peak)

  # merge peaks closer than min_gap, keeping the highest
  keep <- integer(0)
  for (p in peaks[order(-s[peaks])]) {
    if (all(abs(m$t[p] - m$t[keep]) >= min_gap)) keep <- c(keep, p)
  }
  keep <- sort(keep)
  list(peak_t = m$t[keep], peak_share = s[keep], intervals = diff(m$t[keep]), threshold = threshold)
}

# Intervals between consecutive outbreaks in the same cell (months)
rhd_cell_intervals <- function(outbreaks, from = NULL){
  o <- outbreaks
  o$t <- rhd_month_index(o$year, o$month)
  if (!is.null(from)) o <- o[o$t >= from, ]
  o <- o[order(o$x, o$y, o$t), ]
  same_cell <- c(FALSE, o$x[-1] == o$x[-nrow(o)] & o$y[-1] == o$y[-nrow(o)])
  diff(c(NA, o$t))[same_cell]
}

# Outbreak intervals at the chosen scale ("park": between park-wide waves; "cell": within cells)
rhd_intervals <- function(logs, scale = c("park", "cell"), from = NULL, ...){
  scale <- match.arg(scale)
  if (scale == "park") rhd_park_waves(logs$monthly, from = from, ...)$intervals
  else rhd_cell_intervals(logs$outbreaks, from = from)
}

# Distance of the outbreak intervals from the target range, for calibration: 0 when the median is
# inside `target` (months), otherwise how far outside, in units of the range width. Without at
# least two waves the distance is `no_waves` (a large value).
rhd_interval_distance <- function(intervals, target = c(24, 36), no_waves = 5){
  if (length(intervals) == 0) return(no_waves)
  med <- stats::median(intervals)
  width <- diff(target)
  if (med < target[1]) (target[1] - med) / width
  else if (med > target[2]) (med - target[2]) / width
  else 0
}

# Population change after an arrival: mean number of rabbits in the `window` months after the
# arrival month relative to the `window` months before (e.g. -0.6 = a 60% decline). `arrival` is
# "YYYY-MM".
rhd_decline <- function(monthly, arrival, window = 12){
  a <- as.integer(strsplit(arrival, "-")[[1]])
  t0 <- rhd_month_index(a[1], a[2])
  t <- rhd_month_index(monthly$year, monthly$month)
  before <- monthly$rabbits[t < t0 & t >= t0 - window]
  after  <- monthly$rabbits[t > t0 & t <= t0 + window]
  if (length(before) == 0 || length(after) == 0) return(NA_real_)
  mean(after) / mean(before) - 1
}
