#' Choose engine with FFT-friendliness awareness
#'
#' Heuristic that selects \code{"fft"} only when the target linear convolution
#' size is large enough and the implied FFT padding is "friendly" (next 5-smooth
#' sizes do not inflate area too much). Otherwise selects \code{"cpp"}.
#'
#' @param x SpatRaster.
#' @param d Numeric. Window radius/size in map units.
#' @param w Character. Window type.
#' @param fun Character. Summary function.
#' @param min_pixels_fft Numeric. Minimum (n+k-1)*(m+k-1) area to consider FFT.
#' @param min_kernel_fft Numeric. Minimum kernel area to consider FFT.
#' @param max_pad_inflate Numeric. Max allowed area inflation from padding
#'   (e.g., 1.25 allows up to +25%).
#'
#' @return \code{"fft"} or \code{"cpp"}.
#' @keywords internal
#' @importFrom terra nrow ncol
choose_engine_smart <- function(x, d, w = "circle", fun = "mean",
                                min_pixels_fft = 5e6,
                                min_kernel_fft = 81,
                                max_pad_inflate = 1.25) {
  if (!(fun %in% c("mean","sum"))) return("cpp")
  
  nr <- terra::nrow(x); nc <- terra::ncol(x)
  K <- fastfocal_weights(x = x, d = d, w = w, normalize = FALSE)
  kr <- nrow(K); kc <- ncol(K)
  
  tr <- nr + kr - 1; tc <- nc + kc - 1
  next_r <- next_fast_len(tr); next_c <- next_fast_len(tc)
  pad_ratio <- (next_r * next_c) / (tr * tc)
  
  big_enough <- (tr * tc) >= min_pixels_fft || (kr * kc) >= min_kernel_fft
  friendly   <- pad_ratio <= max_pad_inflate
  
  if (big_enough && friendly) "fft" else "cpp"
}
