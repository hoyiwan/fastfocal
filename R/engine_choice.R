#' Choose computation engine (internal)
#'
#' Heuristic to select "fft" vs "cpp" based on raster size, kernel size,
#' and padding cost. Only `mean` and `sum` consider the FFT path. The kernel
#' size is estimated from `d / res(x)` without building the kernel.
#'
#' @param x SpatRaster. Input raster.
#' @param d numeric. Kernel radius/size in map units.
#' @param w character. Window type name (e.g., "circle"); used only to estimate support size.
#' @param fun character. Summary function; FFT considered only for `mean` and `sum`.
#' @param pad character. "auto" pads to next 5-smooth lengths; "none" uses exact sizes.
#' @param min_pixels_fft numeric. Minimum padded conv area to consider FFT (default 5e6).
#' @param min_kernel_fft numeric. Minimum kernel support (cells) to consider FFT (default 81).
#' @param max_pad_inflate numeric. Maximum padding inflation ratio (e.g., 1.25 allows +25%).
#' @return "fft" or "cpp".
#' @keywords internal
#' @noRd

choose_engine_smart <- function(x, d, w = "circle", fun = "mean",
                                pad = "auto",
                                min_pixels_fft = 5e6,
                                min_kernel_fft = 81,
                                max_pad_inflate = 1.25) {
  # Only mean/sum benefit from the FFT path here
  if (!fun %in% c("mean", "sum")) return("cpp")
  
  # Basic geometry
  nr <- as.integer(terra::nrow(x)); if (!is.finite(nr) || nr <= 0) return("cpp")
  nc <- as.integer(terra::ncol(x)); if (!is.finite(nc) || nc <= 0) return("cpp")
  
  resv <- terra::res(x)[1]
  if (!is.finite(resv) || resv <= 0) return("cpp")
  
  # Kernel size in cells (square box that contains the window)
  rad <- max(1L, as.integer(ceiling(d / resv)))
  kr  <- 2L * rad + 1L
  kc  <- kr
  
  # Approximate support size: circle uses pi * r^2, others use full box
  if (is.character(w) && w %in% c("circle", "circular")) {
    k_support <- as.numeric(pi) * (rad ^ 2)
  } else {
    k_support <- kr * kc
  }
  
  # Linear convolution sizes and padded sizes
  want_nr <- nr + kr - 1L
  want_nc <- nc + kc - 1L
  pad_nr  <- if (identical(pad, "auto")) next_fast_len(want_nr) else want_nr
  pad_nc  <- if (identical(pad, "auto")) next_fast_len(want_nc) else want_nc
  
  pad_ratio <- (pad_nr * pad_nc) / (want_nr * want_nc)
  
  # Thresholds for "large enough" and "padding is friendly"
  big_enough <- (want_nr * want_nc) >= min_pixels_fft || k_support >= min_kernel_fft
  friendly   <- pad_ratio <= max_pad_inflate
  
  if (big_enough && friendly) "fft" else "cpp"
}
