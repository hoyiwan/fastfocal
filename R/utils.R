#' Resolve Summary Function for Focal Operations
#'
#' Converts a summary function input (e.g., "mean", "median", "p25") into either
#' a valid string recognized by `terra::focal()` or a custom function defined below.
#'
#' @param fun Character string or function. Summary function to use.
#'
#' @return A character string or function object to be passed to `terra::focal()`.
#' @keywords internal
.makeTextFun <- function(fun) {
  if (is.character(fun)) {
    fun <- tolower(fun)
    
    # Known terra-supported string shortcuts
    if (fun %in% c("mean", "sum", "min", "max", "range", "sd")) {
      return(fun)
    }
    
    # Custom functions mapped by name
    fun_map <- list(
      median = .focal_median,
      mode   = .focal_mode,
      p25    = .focal_p25,
      p75    = .focal_p75
    )
    
    if (fun %in% names(fun_map)) {
      return(fun_map[[fun]])
    }
    
    # Try to retrieve as existing function (last resort)
    if (exists(fun, mode = "function", inherits = TRUE)) {
      return(get(fun, mode = "function", inherits = TRUE))
    } else {
      stop(sprintf("Unknown summary function: '%s'", fun))
    }
  } else if (is.function(fun)) {
    return(fun)
  } else {
    stop("`fun` must be a character string or a function.")
  }
}


#' Focal Median
#' @keywords internal
.focal_median <- function(x, na.rm = TRUE) {
  stats::median(x, na.rm = na.rm)
}

#' Focal Mode
#'
#' Returns the most frequent value. In case of ties, the first one is returned.
#' @keywords internal
.focal_mode <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Focal 25th Percentile
#' @keywords internal
.focal_p25 <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.25, na.rm = na.rm, names = FALSE, type = 7)
}

#' Focal 75th Percentile
#' @keywords internal
.focal_p75 <- function(x, na.rm = TRUE) {
  stats::quantile(x, probs = 0.75, na.rm = na.rm, names = FALSE, type = 7)
}
