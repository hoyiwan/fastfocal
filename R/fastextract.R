#' Fast extraction from raster at point locations or buffer scales
#'
#' @param r SpatRaster. Raster with one or more layers.
#' @param pts SpatVector or sf. Point locations.
#' @param stat Character. Summary statistic (e.g., "mean", "min", "max", "sd", "sum").
#' @param scales Numeric. Buffer radii in meters. Default is 0 (point extraction).
#' @param na.rm Logical. Whether to ignore NA values. Default is TRUE.
#' @param window Character. Window shape: "circular" (default), "rectangular", or "gaussian".
#'
#' @return Matrix of extracted values: rows = points, cols = layers × scales
#' @export
#' @useDynLib fastfocal, .registration = TRUE
#' @importFrom Rcpp evalCpp
fastextract <- function(r, pts, stat = "mean", scales = 0, na.rm = TRUE,
                        window = c("circular", "rectangular", "gaussian")) {
  window <- match.arg(window)
  
  if (inherits(pts, "sf")) pts <- terra::vect(pts)
  stopifnot(inherits(pts, "SpatVector"))
  stopifnot(terra::geomtype(pts) == "points")
  
  layer_names <- names(r)
  if (is.null(layer_names)) {
    layer_names <- paste0("lyr", seq_len(terra::nlyr(r)))
  }
  
  # Fast path: point-only
  if (all(scales == 0)) {
    mats <- lapply(1:terra::nlyr(r), function(i) terra::as.matrix(r[[i]], wide = TRUE))
    coords <- terra::crds(pts)
    ext <- terra::ext(r)
    res <- terra::res(r)
    
    out <- fastextract_multi_cpp(mats, coords,
                                 x_coords = c(ext[1]),
                                 y_coords = c(ext[3]),
                                 res_x = c(res[1]),
                                 res_y = c(res[2]),
                                 na_rm = na.rm)
    colnames(out) <- paste0(layer_names, "_s0")
    return(out)
  }
  
  # Buffer-based extraction path
  vals_list <- list()
  set_fastfocal_progress_handler()
  
  progressr::with_progress({
    p <- progressr::progressor(steps = length(scales))
    
    for (s in scales) {
      mats <- lapply(1:terra::nlyr(r), function(i) terra::as.matrix(r[[i]], wide = TRUE))
      coords <- terra::crds(pts)
      ext <- terra::ext(r)
      res <- terra::res(r)
      
      out <- lapply(seq_along(mats), function(j) {
        fastextract_buffer_cpp(
          mat = mats[[j]],
          coords = coords,
          x_min = ext[1],
          y_min = ext[3],
          res_x = res[1],
          res_y = res[2],
          radius = s,
          stat = stat,
          na_rm = na.rm,
          window = window
        )
      })
      
      out_mat <- do.call(cbind, out)
      colnames(out_mat) <- paste0(layer_names, "_s", s)
      vals_list[[as.character(s)]] <- out_mat
      p()
    }
  })
  
  return(do.call(cbind, unname(vals_list)))
}
