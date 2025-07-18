#' Fast focal smoothing with FFT auto-switching
#'
#' Applies a focal smoothing operation to a SpatRaster using either a C++ (via terra) or FFT backend.
#' Supports a range of window types including rectangle, circle, gaussian, pareto, idw, and exponential.
#'
#' @param x SpatRaster. Input raster layer or multi-layer stack.
#' @param d Numeric. Radius or size of the smoothing window (in map units).
#' @param w Character. Type of window: one of
#'   "rectangle", "circle", "gaussian", "pareto", "idw", "exponential",
#'   "triangular", "cosine", "logistic", "cauchy", "quartic", or "epanechnikov".
#' @param fun Character. Function to apply: one of "mean", "sum", "min", "max", "sd", "median".
#' @param engine Character. Either "auto" (default), "cpp", or "fft" to choose backend.
#' @param na.rm Logical. Remove NAs before applying function.
#' @param ... Additional arguments passed to [terra::focal()] if using the C++ backend.
#'
#' @return A [terra::SpatRaster] of the same dimensions as `x`, containing the focal-smoothed values.
#'
#' @export
#'
#' @importFrom terra rast res values nrow ncol nlyr focal buffer extract vect geomtype
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors
#' @importFrom stats fft median sd
#' 
#' @examples
#' library(terra)
#' x <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000)
#' values(x) <- runif(ncell(x))
#' result <- fastfocal(x, d = 300, w = "circle", fun = "mean")
#' plot(result)

fastfocal <- function(x, d, w = "circle", fun = "mean", engine = "auto", na.rm = TRUE, ...) {
  stopifnot(inherits(x, "SpatRaster"))
  
  # Valid kernel names
  valid_windows <- c("rectangle", "circle", "circular", "gaussian", "Gauss", "pareto", "idw",
                     "exponential", "triangular", "cosine", "logistic", "cauchy", "quartic", "epanechnikov")
  
  # Capture and check unexpected arguments
  dots <- list(...)
  arg_suggestions <- c(
    "window" = "w",
    "windows" = "w",
    "kernel" = "w",
    "kernels" = "w",
    "radius" = "d",
    "function" = "fun",
    "stat" = "fun"
  )
  for (arg in names(dots)) {
    if (arg %in% names(arg_suggestions)) {
      warning(sprintf("Unrecognized argument '%s'. Did you mean '%s'?", arg, arg_suggestions[[arg]]))
    }
  }
  
  # Validate window type
  if (!(w %in% valid_windows)) {
    stop("Invalid window type: '", w, "'. Must be one of: ",
         paste(valid_windows, collapse = ", "))
  }
  
  # Handle multi-layer input
  if (terra::nlyr(x) > 1) {
    out_list <- lapply(1:terra::nlyr(x), function(i) {
      xi <- x[[i]]
      out_i <- fastfocal(xi, d = d, w = w, fun = fun, engine = engine, na.rm = na.rm, ...)
      names(out_i) <- names(x)[i]
      out_i
    })
    return(do.call(c, out_list))
  }
  
  # Auto-select engine
  if (engine == "auto") {
    engine <- choose_engine(x, d, fun = fun)
  }
  
  # Determine whether to normalize the kernel (only for mean)
  normalize_kernel <- (fun %in% c("mean"))
  
  # Build initial kernel (normalized if mean; raw counts if sum)
  kernel <- fastfocal_weights(x = x, d = d, w = w, normalize = normalize_kernel)
  
  if (engine == "fft") {
    mat <- matrix(terra::values(x),
                  nrow = terra::nrow(x),
                  ncol = terra::ncol(x),
                  byrow = TRUE)
    
    # Rebuild the kernel once using the correct normalize flag
    kernel <- fastfocal_weights(x = x, d = d, w = w, normalize = normalize_kernel)
    
    if (fun == "mean") {
      fft_mat <- fft_convolve(mat, kernel, fun = "mean", na.rm = na.rm)
    } else if (fun == "sum") {
      fft_mat <- fft_convolve(mat, kernel, fun = "sum", na.rm = na.rm)
    } else {
      stop("Only 'mean' and 'sum' are supported by the FFT backend. Use engine = 'cpp' for other functions.")
    }
    
    out <- terra::rast(x)
    terra::values(out) <- as.vector(fft_mat)
    names(out) <- paste0("focal_", fun)
    return(out)
    
  } else if (engine == "cpp") {
    fun_eval <- resolve_summary_function(fun)
    return(terra::focal(x, w = kernel, fun = fun_eval, na.rm = na.rm, ...))
    
  } else {
    stop("Unsupported engine: must be 'auto', 'cpp', or 'fft'")
  }
}
