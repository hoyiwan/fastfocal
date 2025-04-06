#' Fast moving window (focal) statistics using circular, rectangular, or Gaussian kernels
#'
#' @param r SpatRaster. Single- or multi-layer raster.
#' @param stat Character. Summary statistic: "mean", "sd", "min", "max", "median", "range", "sum", "p25", "p75".
#' @param radius Numeric. Radius (in map units) for the moving window.
#' @param na.rm Logical. Whether to ignore NA values in the window.
#' @param window Character. Window shape: "circular", "rectangular", "gaussian".
#' @return SpatRaster with focal values for each layer
#' @export
fastfocal <- function(r, stat = "mean", radius = 1, na.rm = TRUE,
                      window = c("circular", "rectangular", "gaussian")) {
  stopifnot(inherits(r, "SpatRaster"))
  window <- match.arg(window)
  
  if (terra::nlyr(r) > 1) {
    cli::cli_warn("Multi-layer input detected. This may take longer to compute.")
  }
  
  ext <- terra::ext(r)
  res <- terra::res(r)
  ncell <- terra::ncell(r)
  
  xy <- terra::xyFromCell(r, 1:ncell)
  pts <- terra::vect(xy, type = "points", crs = terra::crs(r))
  
  layer_names <- names(r)
  if (is.null(layer_names)) {
    layer_names <- paste0("lyr", seq_len(terra::nlyr(r)))
  }
  
  vals <- fastextract(r, pts, stat = stat, scales = radius,
                      na.rm = na.rm, window = window)
  
  # Turn matrix into a SpatRaster with 1+ layers
  layers <- list()
  for (i in seq_len(ncol(vals))) {
    m <- matrix(vals[, i], nrow = terra::nrow(r), ncol = terra::ncol(r), byrow = TRUE)
    rr <- terra::rast(m)
    terra::ext(rr) <- ext
    terra::crs(rr) <- terra::crs(r)
    names(rr) <- colnames(vals)[i]
    layers[[i]] <- rr
  }
  
  out <- do.call(c, layers)
  return(out)
}
