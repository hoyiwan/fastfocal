#' Resolve Summary Function for Focal Operations
#'
#' Converts a summary function input (as a character or function) into either
#' a recognized keyword (e.g., "mean") or an R function for use in `terra::focal()`.
#'
#' Supports standard summary functions including `"mean"`, `"sum"`, `"min"`, `"max"`,
#' `"sd"`, and `"median"`. Custom functions can also be supplied directly.
#'
#' @param fun A character string or function. The name or definition of a summary function.
#'
#' @return A character string or function that can be passed to [terra::focal()].
#' @export
#'
#' @importFrom stats median sd
resolve_summary_function <- function(fun) {
  valid_strings <- c("mean", "sum", "min", "max", "sd", "median")
  
  if (is.character(fun)) {
    fun <- tolower(fun)
    
    if (fun %in% valid_strings) {
      return(fun)
    }
    
    if (exists(fun, mode = "function", inherits = TRUE)) {
      return(get(fun, mode = "function", inherits = TRUE))
    }
    
    stop(sprintf("Unsupported summary function: '%s'. Must be one of: %s",
                 fun, paste(valid_strings, collapse = ", ")))
    
  } else if (is.function(fun)) {
    return(fun)
  } else {
    stop("`fun` must be a character string or a function.")
  }
}
