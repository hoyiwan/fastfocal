#' Resolve summary function (internal)
#'
#' Converts a user-supplied summary selector into a value suitable for
#' `terra::focal()`: either one of the keywords "mean","sum","min","max","sd","median",
#' or an actual R function object.
#'
#' @param fun character or function. A single function name or a function.
#' @return character (one of the six keywords) or a function object.
#' @keywords internal
#' @noRd
resolve_summary_function <- function(fun) {
  # If a function was supplied, accept it directly.
  if (is.function(fun)) return(fun)
  
  # Otherwise require a single string.
  if (!is.character(fun) || length(fun) != 1L) {
    stop("`fun` must be a single character string or a function.")
  }
  
  key <- tolower(trimws(fun))
  allowed <- c("mean", "sum", "min", "max", "sd", "median")
  
  # Prefer the well-known keywords — terra accepts these strings.
  if (key %in% allowed) {
    return(key)
  }
  
  # Fallback: if a function of that name exists in scope, return it.
  if (exists(key, mode = "function", inherits = TRUE)) {
    return(get(key, mode = "function", inherits = TRUE))
  }
  
  stop(
    "Unsupported summary function: '", fun,
    "'. Use one of: ", paste(allowed, collapse = ", "),
    "; or supply a function."
  )
}
