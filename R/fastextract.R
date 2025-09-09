#' Fast raster extraction at points (buffered)
#'
#' Extracts summary statistics from a SpatRaster at point locations,
#' optionally using buffered extraction with custom kernel windows.
#'
#' - If `d > 0`, a buffer of radius `d` (map units) is created around each point
#'   and the summary is computed over raster cells intersecting the buffer.
#' - If `d == 0`, values are taken at the point locations (no buffering).
#' - If `y` is a polygon layer, the summary is computed over polygon areas.
#'
#' @param x SpatRaster. Input raster (single- or multi-layer).
#' @param y SpatVector. Points or polygons.
#' @param d numeric or numeric vector. Buffer radius/radii in map units.
#' @param w character. Window type for the buffer kernel when `d > 0`
#'   (currently passed through to \pkg{terra}; e.g., "circle", "rectangle").
#' @param fun character or function. Summary function: "mean", "sum", "min",
#'   "max", "sd", or "median"; or a user function.
#' @param na.rm logical. Whether to remove NAs when computing summaries.
#'
#' @return A data.frame of extracted values. When `d` has multiple values,
#'   rows are stacked by scale with a `scale_m` column indicating the radius.
#'
#' @export
#'
#' @importFrom terra buffer extract geomtype rast res values nrow ncol nlyr vect
#' @importFrom stats median sd
#'
#' @examples
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
#' terra::values(r) <- seq_len(terra::ncell(r))
#'
#' pts <- terra::vect(
#'   matrix(c(10, 10,
#'            50, 50), ncol = 2, byrow = TRUE),
#'   type = "points",
#'   crs  = terra::crs(r)
#' )
#'
#' # Mean over a 20-unit circular neighborhood around each point
#' res <- fastextract(r, pts, d = 20, w = "circle", fun = "mean")
#' head(res)

fastextract <- function(x, y, d = 0, w = "circle", fun = "mean", na.rm = TRUE) {
  if (!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
  if (!inherits(y, "SpatVector")) stop("y must be a SpatVector")
  
  fun_eval <- switch(
    if (is.character(fun)) tolower(fun) else "",
    mean   = base::mean,
    sum    = base::sum,
    min    = base::min,
    max    = base::max,
    sd     = stats::sd,
    median = stats::median,
    { if (is.function(fun)) fun else NULL }
  )
  
  do_one <- function(d_single) {
    if (terra::geomtype(y)[1] == "polygons") {
      if (is.null(fun_eval)) stop("Unsupported summary function for polygons: ", fun)
      df <- terra::extract(x, y, fun = fun_eval, na.rm = na.rm, ID = FALSE)
      
    } else if (!is.null(d_single) && d_single > 0) {
      if (is.null(fun_eval)) stop("Unsupported summary function for buffered extraction: ", fun)
      buff <- terra::buffer(y, width = d_single)
      df <- terra::extract(x, buff, fun = fun_eval, na.rm = na.rm, ID = FALSE)
      
    } else {
      vals <- terra::extract(x, y, ID = FALSE)
      if (is.null(fun_eval)) {
        df <- as.data.frame(t(apply(vals, 1, function(row) {
          vapply(row, function(z) fun(z, na.rm = na.rm), numeric(1))
        })))
        names(df) <- names(x)
      } else {
        df <- as.data.frame(vals)
      }
    }
    df <- as.data.frame(df)
    if ("ID" %in% names(df)) df <- df[, setdiff(names(df), "ID"), drop = FALSE]
    df
  }
  
  if (length(d) > 1) {
    out_list <- lapply(d, function(di) {
      df <- do_one(di)
      df <- as.data.frame(df)
      # prepend scale_m without passing check.names
      df <- cbind(scale_m = di, df)
      df
    })
    return(do.call(rbind, out_list))
  }
  
  do_one(d)
}
