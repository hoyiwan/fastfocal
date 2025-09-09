#' Fast focal smoothing with FFT auto-switching
#'
#' Applies a focal operation to a `SpatRaster` using either a 'C++' backend
#' (via \pkg{terra}) or an 'FFT' backend. Window types include rectangle, circle,
#' gaussian, pareto, idw, exponential, triangular, cosine, logistic, cauchy,
#' quartic, epanechnikov, or you may pass a numeric matrix as the kernel.
#'
#' The 'FFT' backend uses masked convolution with proper NA semantics and can
#' pad to "5-smooth" sizes for stable speed. With `engine = "auto"`, the function
#' chooses between 'C++' and 'FFT' based on a simple window-size heuristic.
#'
#' @param x SpatRaster. Input raster (1+ layers).
#' @param d numeric. Radius/size in map units (ignored if `w` is a matrix).
#' @param w character or numeric matrix. Window type, or a custom kernel matrix.
#' @param fun character. One of "mean","sum","min","max","sd","median".
#' @param engine character. "auto" (default), "cpp", or "fft".
#' @param na.rm logical. Remove NAs before applying the summary function.
#' @param na.policy character. "omit" (default) leaves NA centers as NA; "all" fills
#'   centers when neighbors exist (FFT path respects this; C++ path emulates center
#'   handling after the call).
#' @param pad character. "none" or "auto" (pad to next 5-smooth sizes for FFT).
#' @param ... Extra args to [terra::focal()] for the 'C++' path.
#'
#' @return [terra::SpatRaster] with the same geometry as `x`.
#' @export
#' @importFrom terra rast res values nrow ncol nlyr focal ncell crs ext
#' @importFrom stats median sd
#'
#' @examples
#' set.seed(1)
#' r <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 12)
#' terra::values(r) <- stats::runif(terra::ncell(r))
#'
#' # Mean with a small circular window (d is in map units; here res = 1)
#' m_circ <- fastfocal(r, d = 2, w = "circle", fun = "mean")
#'
#' # Same idea using a custom 3x3 box kernel (uniform mean)
#' k3 <- matrix(1, 3, 3)
#' m_box <- fastfocal(r, w = k3, fun = "mean")
#'
#' # Tiny numeric summaries (keeps examples fast & quiet for CRAN)
#' as.numeric(terra::global(m_circ, "mean", na.rm = TRUE))
#' as.numeric(terra::global(m_box,  "mean", na.rm = TRUE))
fastfocal <- function(
    x, d, w = "circle", fun = "mean",
    engine = "auto", na.rm = TRUE,
    na.policy = c("omit","all"),
    pad = c("none","auto"),
    ...
) {
  stopifnot(inherits(x, "SpatRaster"))
  na.policy <- match.arg(na.policy)
  pad       <- match.arg(pad)
  
  # Valid named windows (implemented in fastfocal_weights)
  is_kernel_matrix <- is.matrix(w)
  valid_windows <- c(
    "rectangle","circle","circular","gaussian","Gauss","pareto","idw",
    "exponential","triangular","cosine","logistic","cauchy","quartic","epanechnikov"
  )
  if (!is_kernel_matrix && !(w %in% valid_windows)) {
    stop("Invalid window type: '", w, "'. Must be one of: ",
         paste(valid_windows, collapse = ", "),
         " , or pass a numeric matrix as `w`.")
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
  
  # --- choose engine (force C++ for non-mean/sum) ---
  if (identical(engine, "auto")) {
    if (!fun %in% c("mean","sum")) {
      engine <- "cpp"
    } else {
      engine <- .ff_choose_engine_smart(x, d, w = w, fun = fun)
    }
  }
  if (identical(engine, "fft") && !fun %in% c("mean","sum")) {
    # user forced FFT but it's unsupported for this reducer → fall back silently
    engine <- "cpp"
  }
  
  # --- resolve kernel as a numeric matrix (UN-normalized) ---
  kernel_mat <- .ff_kernel_matrix(x, d, w)
  
  # -------------------- FFT path --------------------
  if (identical(engine, "fft")) {
    # Raster -> matrix in terra orientation
    xm <- as.matrix(x, wide = TRUE)
    
    # FFT core: fft_convolve_masked() returns t(out)
    outm_raw <- fft_convolve_masked(
      x         = xm,
      kernel    = kernel_mat,
      fun       = fun,
      na.rm     = na.rm,
      na.policy = na.policy,
      pad       = pad
    )
    outm <- t(outm_raw)  # bring back to as.matrix(..., wide=TRUE)
    
    # Matrix -> raster with exact cell order
    out <- .ff_mat_to_rast_exact(outm, x)
    names(out) <- paste0("focal_", fun)
    return(out)
  }
  
  # -------------------- C++ (terra) path --------------------
  fun <- match.arg(fun, c("mean","sum","min","max","sd","median"))
  fun_eval <- .ff_resolve_summary_function(fun)
  
  # Terra expectations:
  # - mean: weights normalized to sum 1
  # - sum : weights as-is
  # - others: 0/1 footprint
  if (fun %in% c("mean","sum")) {
    ksum <- sum(kernel_mat)
    if (identical(fun, "mean")) {
      if (ksum == 0) stop("Kernel sum is zero; cannot normalize for 'mean'.")
      kernel_cpp <- kernel_mat / ksum
    } else {
      kernel_cpp <- kernel_mat
    }
  } else {
    kernel_cpp <- (kernel_mat != 0) * 1
  }
  
  out <- terra::focal(x, w = kernel_cpp, fun = fun_eval, na.rm = na.rm, ...)
  
  # Emulate na.policy center handling (terra always leaves NA when center is NA).
  if (identical(na.policy, "omit")) {
    center_na <- is.na(terra::values(x))
    v <- terra::values(out); v[center_na] <- NA_real_; terra::values(out) <- v
  }
  out
}

# ----- helpers (internal) -----

.ff_resolve_summary_function <- function(fun) {
  switch(fun,
         mean   = base::mean,
         sum    = base::sum,
         min    = base::min,
         max    = base::max,
         sd     = stats::sd,
         median = stats::median)
}

.ff_choose_engine_smart <- function(x, d, w, fun) {
  # Non-mean/sum always use C++.
  if (!fun %in% c("mean","sum")) return("cpp")
  
  # simple heuristic: if footprint >= 15x15 -> FFT, else C++
  if (is.matrix(w)) {
    kr <- nrow(w); kc <- ncol(w)
  } else {
    K  <- fastfocal_weights(x, d, w, normalize = FALSE, plot = FALSE)
    kr <- nrow(K); kc <- ncol(K)
  }
  if ((kr * kc) >= 225L) "fft" else "cpp"
}

.ff_kernel_matrix <- function(x, d, w) {
  if (is.matrix(w)) {
    K <- w
  } else if (is.numeric(w) && length(dim(w)) == 2L) {
    K <- w
  } else if (is.character(w)) {
    if (missing(d) || is.null(d))
      stop("When w is a character window name (e.g. 'circle','gaussian'), you must supply d.")
    # >>> IMPORTANT: normalize=TRUE to mirror terra::focalMat() for named kernels
    K <- fastfocal_weights(x = x, d = d, w = w, normalize = TRUE, plot = FALSE)
  } else {
    stop("`w` must be a numeric matrix or a supported window name (e.g. 'circle','gaussian').")
  }
  K <- as.matrix(K); storage.mode(K) <- "double"; K
}

# exact cell-order mapping back to raster
.ff_mat_to_rast_exact <- function(M, like) {
  stopifnot(is.matrix(M))
  out <- terra::rast(nrows = nrow(like), ncols = ncol(like),
                     ext = terra::ext(like), crs = terra::crs(like))
  idx <- terra::rast(nrows = nrow(like), ncols = ncol(like),
                     ext = terra::ext(like), crs = terra::crs(like),
                     vals = 1:terra::ncell(like))
  idx_m <- as.matrix(idx, wide = TRUE)
  vals <- rep(NA_real_, terra::ncell(like))
  vals[as.vector(idx_m)] <- as.vector(M)
  terra::values(out) <- vals
  out
}

# alias to avoid renaming
mat_to_rast_exact <- .ff_mat_to_rast_exact
