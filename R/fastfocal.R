#' Fast focal smoothing with FFT auto-switching
#'
#' Applies a focal operation to a `SpatRaster` using either a 'C++' backend
#' (via \pkg{terra}) or an 'FFT' backend. Window types include rectangle, circle,
#' gaussian, pareto, idw, exponential, triangular, cosine, logistic, cauchy,
#' quartic, epanechnikov - or you may pass a numeric matrix as the kernel.
#'
#' The 'FFT' backend uses masked convolution with proper NA semantics and can
#' pad to "5-smooth" sizes for stable speed. With `engine = "auto"`, the function
#' chooses between 'C++' and 'FFT' based on scale and FFT-friendliness.
#'
#' @param x SpatRaster. Input raster (1+ layers).
#' @param d numeric. Radius/size in map units (ignored if `w` is a matrix).
#' @param w character or numeric matrix. Window type, or a custom kernel matrix.
#' @param fun character. One of "mean","sum","min","max","sd","median".
#' @param engine character. "auto" (default), "cpp", or "fft".
#' @param na.rm logical. Remove NAs before applying the summary function.
#' @param na.policy character. "omit" (default) leaves NA centers as NA; "all" fills centers when neighbors exist.
#' @param pad character. "none" or "auto" (pad to next 5-smooth sizes for FFT).
#' @param ... Extra args to [terra::focal()] for the 'C++' path.
#'
#' @return [terra::SpatRaster] with the same geometry as `x`.
#' @export
#' @importFrom terra rast res values nrow ncol nlyr focal buffer extract vect geomtype ncell
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors
#' @importFrom stats fft median sd
#'
#' @examples
#' # Tiny, fast example (no plotting, no side effects)
#' library(terra)
#' r <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
#' values(r) <- runif(ncell(r))
#' out_mean <- fastfocal(r, d = 15, w = "circle", fun = "mean")
#' out_sum  <- fastfocal(r, d = 15, w = "circle", fun = "sum")
fastfocal <- function(
    x, d, w = "circle", fun = "mean",
    engine = "auto", na.rm = TRUE,
    na.policy = c("omit","all"),
    pad = c("none","auto"),
    ...
) {
  stopifnot(inherits(x, "SpatRaster"))
  na.policy <- match.arg(na.policy)
  pad <- match.arg(pad)
  
  # Allow a custom matrix kernel
  is_kernel_matrix <- is.matrix(w)
  
  # Valid named windows
  valid_windows <- c(
    "rectangle","circle","circular","gaussian","Gauss","pareto","idw",
    "exponential","triangular","cosine","logistic","cauchy","quartic","epanechnikov"
  )
  
  # Suggest common misnamed args
  dots <- list(...)
  arg_suggestions <- c(
    "window"="w","windows"="w","kernel"="w","kernels"="w",
    "radius"="d","function"="fun","stat"="fun"
  )
  for (arg in names(dots)) {
    if (arg %in% names(arg_suggestions)) {
      warning(sprintf("Unrecognized argument '%s'. Did you mean '%s'?", arg, arg_suggestions[[arg]]))
    }
  }
  
  if (!is_kernel_matrix && !(w %in% valid_windows)) {
    stop("Invalid window type: '", w, "'. Must be one of: ",
         paste(valid_windows, collapse = ", "),
         " - or pass a numeric matrix as `w`.")
  }
  
  # Multi-layer: process layers independently
  if (terra::nlyr(x) > 1) {
    out_list <- lapply(seq_len(terra::nlyr(x)), function(i) {
      xi <- x[[i]]
      out_i <- fastfocal(
        xi, d = d, w = w, fun = fun,
        engine = engine, na.rm = na.rm,
        na.policy = na.policy, pad = pad, ...
      )
      names(out_i) <- names(x)[i]
      out_i
    })
    return(do.call(c, out_list))
  }
  
  # Engine auto-choice
  if (identical(engine, "auto")) {
    engine <- choose_engine_smart(x, d, w = if (is_kernel_matrix) "rectangle" else w, fun = fun)
  }
  
  # Build kernel (unnormalized) - use matrix as given
  if (is_kernel_matrix) {
    kernel <- w
    if (!is.numeric(kernel)) stop("Custom kernel matrix must be numeric.")
    if (any(is.na(kernel))) stop("Custom kernel matrix must not contain NA.")
  } else {
    kernel <- fastfocal_weights(x = x, d = d, w = w, normalize = FALSE)
  }
  ksum <- sum(kernel)
  
  if (identical(engine, "fft")) {
    if (!fun %in% c("mean","sum")) stop("FFT backend supports only 'mean' and 'sum'.")
    mat <- matrix(terra::values(x), nrow = terra::nrow(x), ncol = terra::ncol(x), byrow = TRUE)
    fft_mat <- fft_convolve_masked(
      mat, kernel, fun = fun, na.rm = na.rm,
      na.policy = na.policy, pad = pad
    )
    out <- terra::rast(x)
    terra::values(out) <- as.vector(fft_mat)
    names(out) <- paste0("focal_", fun)
    return(out)
    
  } else if (identical(engine, "cpp")) {
    fun_eval <- resolve_summary_function(fun)
    kernel_cpp <- if (identical(fun, "mean")) kernel / ksum else kernel
    out <- terra::focal(x, w = kernel_cpp, fun = fun_eval, na.rm = na.rm, ...)
    if (identical(na.policy, "omit")) {
      center_na <- is.na(terra::values(x))
      v <- terra::values(out); v[center_na] <- NA_real_; terra::values(out) <- v
    }
    return(out)
    
  } else {
    stop("Unsupported engine: must be 'auto', 'cpp', or 'fft'")
  }
}
