#' Select the computation engine ("cpp" vs "fft") for fastfocal()
#'
#' Heuristic (pixel-only, resolution-independent):
#' - Only additive reducers (\code{fun} in \code{c("mean","sum")}) are eligible for FFT.
#' - Let \eqn{r_\mathrm{cells}} be the kernel radius in cells (derived from \code{w} or \code{d/res(x)}).
#'   If \eqn{r_\mathrm{cells} \le 14}, use \code{"cpp"}.
#' - Otherwise (\eqn{r_\mathrm{cells} \ge 15}), compute the padded linear-convolution grid:
#'   \deqn{\textit{want\_nr} = n_\mathrm{row} + (2r_\mathrm{cells}+1) - 1,\quad
#'         \textit{want\_nc} = n_\mathrm{col} + (2r_\mathrm{cells}+1) - 1}
#'   Then set \eqn{\textit{pad\_nr}}, \eqn{\textit{pad\_nc}} to the next 5-smooth lengths
#'   when \code{pad == "auto"} (via \code{next_fast_len()}), otherwise leave as-is.
#'   Switch to \code{"fft"} only if
#'   \deqn{\textit{pad\_area} = \textit{pad\_nr}\times\textit{pad\_nc} \;>\; \max(16{,}000,\; 0.6\, n_\mathrm{row} n_\mathrm{col}).}
#'
#' Equality at thresholds stays on \code{"cpp"}.
#'
#' @param x   A \code{terra::SpatRaster}.
#' @param d   Numeric. Window radius/size in map units (used only to derive cells when \code{w} is named).
#' @param w   Character window name (e.g., \code{"circle"}) or a numeric kernel matrix; used to get the discrete footprint.
#' @param fun Character. Summary function. Only \code{"mean"} and \code{"sum"} may use FFT.
#' @param pad Character. \code{"auto"} pads each dimension to a 5-smooth length; \code{"none"} uses exact sizes.
#'
#' @return A length-1 character: \code{"cpp"} or \code{"fft"}.
#' @keywords internal
#' @noRd
choose_engine_smart <- function(x, d, w = "circle", fun = "mean", pad = "auto") {
  # 1) Reducer gate --------------------------------------------------------
  if (!fun %in% c("mean","sum")) return("cpp")
  
  # 2) Validity guard ------------------------------------------------------
  nr <- as.integer(terra::nrow(x)); nc <- as.integer(terra::ncol(x))
  if (!is.finite(nr) || !is.finite(nc) || nr <= 0L || nc <= 0L) return("cpp")
  resx <- terra::res(x)[1]
  if (!is.finite(resx) || resx <= 0) return("cpp")
  
  # 3) Discrete kernel geometry (respect aliases; NEVER normalize here) ---
  if (is.matrix(w) || (is.numeric(w) && length(dim(w)) == 2L)) {
    kr <- as.integer(nrow(w)); kc <- as.integer(ncol(w))
  } else {
    K  <- fastfocal_weights(x, d, w, normalize = FALSE, plot = FALSE)
    kr <- as.integer(nrow(K)); kc <- as.integer(ncol(K))
  }
  if (!is.finite(kr) || !is.finite(kc) || kr <= 0L || kc <= 0L) return("cpp")
  r_cells <- (max(kr, kc) - 1L) %/% 2L
  
  # 4) Hard pixel-radius switch -------------------------------------------
  if (r_cells <= 14L) return("cpp")  # ≤14 px radius -> C++
  
  # 5) Padded convolution grid area ---------------------------------------
  kw <- 2L * r_cells + 1L
  want_nr <- nr + kw - 1L
  want_nc <- nc + kw - 1L
  pad_nr  <- if (identical(pad, "auto")) next_fast_len(want_nr) else want_nr
  pad_nc  <- if (identical(pad, "auto")) next_fast_len(want_nc) else want_nc
  pad_area <- as.numeric(pad_nr) * as.numeric(pad_nc)
  
  # 6) Adaptive size gate --------------------------------------------------
  ncell <- as.numeric(nr) * as.numeric(nc)
  T <- max(16000.0, 0.6 * ncell)  # small absolute floor + relative term
  
  if (!is.finite(pad_area) || pad_area <= T) return("cpp")
  
  # 7) Otherwise, FFT ------------------------------------------------------
  "fft"
}
