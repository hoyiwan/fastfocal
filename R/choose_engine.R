#' Choose Engine for Fast Focal Operation
#'
#' Selects "fft" or "cpp" based on raster size, kernel size, and summary function.
#'
#' @param r SpatRaster. The input raster.
#' @param d Numeric. Kernel distance/radius.
#' @param fun Character. Summary function (e.g., "mean", "sum", "sd").
#' @param threshold Numeric. Heuristic threshold for fft (default: 1e6).
#'
#' @return Character string: "fft" or "cpp"
#' 
#' @importFrom terra ncell
#' 
#' @keywords internal
choose_engine <- function(r, d, fun = "mean", threshold = 1e6) {
  # Only mean and sum are supported by FFT
  if (!(fun %in% c("mean", "sum"))) {
    return("cpp")
  }
  
  ncells <- terra::ncell(r)
  if ((d^2) * ncells >= threshold) {
    return("fft")
  } else {
    return("cpp")
  }
}
