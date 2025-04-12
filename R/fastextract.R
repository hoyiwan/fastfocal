#' Fast raster extraction at points (buffered)
#'
#' Extracts summary statistics from a SpatRaster at point locations,
#' optionally using buffered extraction with custom kernel windows.
#'
#' @param x SpatRaster. Input raster.
#' @param y SpatVector. Points, polygons, etc.
#' @param d Numeric or vector. Buffer radius/radii in map units.
#' @param w Character. Window type for kernel: "circle", "rectangle", etc.
#' @param fun Character or function. Summary function: "mean", "sum", "min", "max", "sd", "median".
#' @param na.rm Logical. Whether to remove NAs.
#'
#' @return A data.frame of extracted values, one row per point per scale.
#' @export
fastextract <- function(x, y, d = 0, w = "circle", fun = "mean", na.rm = TRUE) {
  if (!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
  if (!inherits(y, "SpatVector")) stop("y must be a SpatVector")
  
  if (length(d) > 1) {
    out_list <- lapply(d, function(di) {
      df <- fastextract(x, y, d = di, w = w, fun = fun, na.rm = na.rm)
      if ("ID" %in% names(df)) {
        df <- df[, !(names(df) %in% "ID"), drop = FALSE]
      }
      df$scale_m <- di
      df <- df[, c("scale_m", setdiff(names(df), "scale_m"))]
      return(df)
    })
    return(do.call(rbind, out_list))
  }
  if (geomtype(y)[1] == "polygons") {
    # Special case: user supplied a polygon layer
    df <- terra::extract(x, y, fun = match.fun(fun), na.rm = na.rm, ID = TRUE)
    
    } else if (!is.null(d) && d > 0) {
      buff <- terra::buffer(y, width = d)
      df <- terra::extract(x, buff, fun = match.fun(fun), na.rm = na.rm)
    
      } else {
        df <- terra::extract(x, y)[, -1, drop = FALSE]  # drop ID
        fun_eval <- switch(
          tolower(fun),
          mean   = base::mean,
          sum    = base::sum,
          min    = base::min,
          max    = base::max,
          sd     = stats::sd,
          median = stats::median,
          stop("Unsupported summary function: ", fun)
          )
    
    # Apply summary function per row, per column
    df <- as.data.frame(t(apply(df, 1, function(row) {
      vapply(row, function(x) fun_eval(x, na.rm = na.rm), numeric(1))
    })))
    names(df) <- names(x)
  }
  
  return(df)
}