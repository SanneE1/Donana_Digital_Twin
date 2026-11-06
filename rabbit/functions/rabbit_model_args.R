# Command line of the rabbit model (Pascal program, settings interface 2).
#
# Every setting is passed as --name=value (see Pascal_program/model_settings.md):
#   code defaults < rates file (--rates_file) < --settings files < --name=value
# Unknown names stop the program. The same helper is used by model_fn(),
# run_simulation() and the forecasting component, so they all build the call
# the same way.

RABBIT_SETTINGS_INTERFACE <- 2

# Named list -> c("--name=value", ...). NULL entries are left out; numbers are
# written with "." and without scientific notation; logicals as 0/1. Paths
# are quoted, so they may contain spaces.
rabbit_model_args <- function(settings){
  settings <- settings[!vapply(settings, is.null, logical(1))]
  if (length(settings) == 0) return(character(0))
  if (is.null(names(settings)) || any(!nzchar(names(settings)))) stop("settings must be a named list")
  if (anyDuplicated(names(settings))) {
    stop("settings given twice: ", paste(unique(names(settings)[duplicated(names(settings))]), collapse = ", "))
  }
  if (any(lengths(settings) != 1)) stop("every setting must be a single value")
  values <- vapply(settings, function(v){
    if (is.logical(v)) v <- as.integer(v)
    if (is.numeric(v)) format(v, scientific = FALSE, trim = TRUE, digits = 15) else as.character(v)
  }, character(1))
  shQuote(paste0("--", names(settings), "=", values), type = if (.Platform$OS.type == "windows") "cmd" else "sh")
}

# Calibrated parameters -> model settings, for the parameters that the Pascal
# program reads (not the habitat parameters, which go into the K maps):
#   R_lambda, R_sigma     -> lambda, sigma (dens_opt keeps its name, rounded)
#   log10_<name>          -> <name> = 10^value (e.g. log10_RHD_spark)
#   RHDMort               -> RHDMort_L = RHDMort_U = value (schedule mode, one mortality)
#   RHD_..., RHDFreq_...  -> as they are
# `p` is a named list or a one-row data frame; other names are ignored.
model_param_settings <- function(p){
  p <- as.list(p)
  out <- list()
  if (!is.null(p$dens_opt)) out$dens_opt <- round(p$dens_opt)
  if (!is.null(p$R_lambda)) out$lambda <- p$R_lambda
  if (!is.null(p$R_sigma))  out$sigma  <- p$R_sigma
  for (nm in grep("^(log10_)?RHD", names(p), value = TRUE)) {
    value <- p[[nm]]
    if (startsWith(nm, "log10_")) { nm <- sub("^log10_", "", nm); value <- 10^value }
    if (nm == "RHDMort") {
      out$RHDMort_L <- value
      out$RHDMort_U <- value
    } else {
      out[[nm]] <- value
    }
  }
  out
}

# Stops with a clear message when the executable is older than the named
# argument interface. Call it once per session, not per run.
check_rabbit_model <- function(exe){
  out <- tryCatch(suppressWarnings(system2(exe, "--version", stdout = TRUE, stderr = TRUE)),
                  error = function(e) conditionMessage(e))
  ok <- any(grepl(paste0("^settings_interface ", RABBIT_SETTINGS_INTERFACE, "$"), trimws(out)))
  if (!ok) {
    stop("The rabbit model '", exe, "' does not support named arguments (settings interface ",
         RABBIT_SETTINGS_INTERFACE, "). Use a build of Pascal_program from branch ",
         "named-model-settings or later. Output of --version: ",
         paste(utils::head(out, 3), collapse = " | "), call. = FALSE)
  }
  invisible(TRUE)
}
