#' 2D FFT convolution (compatibility wrapper)
#'
#' Convenience wrapper that performs 2D convolution via FFT using the same
#' NA semantics and padding strategy as `fastfocal(engine = "fft")`.
#'
#' Internally this delegates to the package's masked FFT implementation
#' (box gating when `na.rm = FALSE`, support-based averaging when
#' `na.rm = TRUE`), and pads to next 5-smooth lengths for stable speed.
#'
#' @param x numeric matrix. Input values (rows x cols).
#' @param kernel numeric matrix. Unnormalized kernel (e.g., 0/1 circle).
#' @param fun character. "mean" or "sum".
#' @param na.rm logical. If `TRUE`, ignore NAs under the kernel; if `FALSE`,
#'   require a full valid box (terra-compatible gating).
#' @return numeric matrix (rows x cols), aligned to `x` and oriented like
#'   `terra::values`, i.e., ready to be assigned back after `as.vector()`.
#' @export
#' @examples
#' m <- matrix(1:9, 3, 3)
#' k <- matrix(1, 3, 3)
#' fft_convolve(m, k, fun = "sum")
#'
#' # Mean with missing center:
#' m_na <- m; m_na[2, 2] <- NA
#' fft_convolve(m_na, k, fun = "mean", na.rm = TRUE)
fft_convolve <- function(x, kernel, fun = "mean", na.rm = TRUE) {
  # Basic checks and normalization of inputs
  x <- as.matrix(x)
  kernel <- as.matrix(kernel)
  
  if (!fun %in% c("mean", "sum")) {
    stop("`fun` must be either 'mean' or 'sum'.")
  }
  
  # Delegate to the masked FFT with stable padding and terra-compatible NA policy.
  # We use na.policy = "omit" here to mirror typical terra usage and the package defaults.
  fft_convolve_masked(
    x       = x,
    kernel  = kernel,
    fun     = fun,
    na.rm   = na.rm,
    na.policy = "omit",
    pad     = "auto"
  )
}
