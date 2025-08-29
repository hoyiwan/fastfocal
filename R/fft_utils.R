#' Masked FFT convolution with next-fast-length padding
#'
#' Performs 2D convolution via FFT with correct NA semantics by using a
#' value mask and (optionally) padding each dimension to a "5-smooth"
#' length (product of 2, 3, 5) for stable performance. Results are cropped
#' to the original raster size and oriented to match \pkg{terra} values.
#'
#' @param mat Numeric matrix. Input raster values (rows × cols).
#' @param kernel Numeric matrix. Focal weights (UNnormalized). Any non-binary
#'   kernel is supported; for binary footprints (e.g., circle/rectangle) this
#'   is typically 0/1.
#' @param fun Character. One of \code{"mean"} or \code{"sum"}.
#' @param na.rm Logical. If \code{TRUE}, ignore NAs; if \code{FALSE}, require
#'   a fully valid window (any NA in window yields NA).
#' @param na.policy Character. \code{"omit"} leaves NA centers as NA;
#'   \code{"all"} fills NA centers when neighbors exist (gap-filling).
#' @param pad Character. \code{"none"} or \code{"auto"} (pad to next 5-smooth
#'   sizes). \code{"auto"} stabilizes speed on awkward sizes.
#'
#' @return Numeric matrix (rows × cols) aligned with \code{mat}, transposed
#'   at the end to match \pkg{terra} value ordering used by \code{fastfocal()}.
#'
#' @keywords internal
#' @importFrom stats fft
fft_convolve_masked <- function(mat, kernel, fun = c("mean","sum"),
                                na.rm = TRUE, na.policy = c("omit","all"),
                                pad = c("none","auto")) {
  fun <- match.arg(fun)
  na.policy <- match.arg(na.policy)
  pad <- match.arg(pad)
  
  nr <- nrow(mat); nc <- ncol(mat)
  kr <- nrow(kernel); kc <- ncol(kernel)
  
  # Value matrix (NA->0) and 0/1 validity mask for input
  center_valid <- is.finite(mat)
  V <- mat; V[!center_valid] <- 0
  M <- matrix(0, nr, nc); M[center_valid] <- 1
  
  # Linear convolution target and padding
  tr <- nr + kr - 1; tc <- nc + kc - 1
  if (pad == "auto") { pr <- next_fast_len(tr); pc <- next_fast_len(tc) } else { pr <- tr; pc <- tc }
  
  pad_to <- function(A, r, c) { B <- matrix(0, r, c); B[seq_len(nrow(A)), seq_len(ncol(A))] <- A; B }
  Vp <- pad_to(V, pr, pc)
  Mp <- pad_to(M, pr, pc)
  Kp <- pad_to(kernel, pr, pc)
  
  # --- Two kernels: weighted for the sum, and NA-scope for the gating ---
  # weighted kernel (actual operation)
  G <- fft(Kp)
  
  # NA scope kernel:
  # - terra behavior: when na.rm = FALSE, ANY NA in the **whole box** (kr x kc) invalidates
  # - when na.rm = TRUE, averages/sums ignore NAs only over the **support** of kernel (w != 0)
  if (na.rm) {
    K_scope <- (kernel != 0) * 1  # support mask
    scope_mass <- sum(K_scope)
  } else {
    K_scope <- matrix(1, kr, kc)  # FULL BOX
    scope_mass <- kr * kc
  }
  Ksp <- pad_to(K_scope, pr, pc)
  Gscope <- fft(Ksp)
  
  # FFTs
  F <- fft(Vp)
  H <- fft(Mp)
  
  # Convolutions
  S_full <- Re(fft(F * G,      inverse = TRUE)) / (pr * pc)  # weighted sum over available values
  C_full <- Re(fft(H * Gscope, inverse = TRUE)) / (pr * pc)  # count of available cells in NA-scope
  
  # Crop (linear conv) and center to original size
  S_lin <- S_full[seq_len(tr), seq_len(tc)]
  C_lin <- C_full[seq_len(tr), seq_len(tc)]
  r_off <- floor((kr - 1) / 2); c_off <- floor((kc - 1) / 2)
  r_idx <- (1:nr) + r_off; c_idx <- (1:nc) + c_off
  S <- S_lin[r_idx, c_idx, drop = FALSE]
  C <- C_lin[r_idx, c_idx, drop = FALSE]
  
  # Tolerance and clamping for FFT roundoff
  tol <- 1e-8 * max(1, scope_mass)
  Cc <- pmin(pmax(C, 0), scope_mass)
  
  out <- matrix(NA_real_, nr, nc)
  
  if (fun == "sum") {
    if (na.rm) {
      # sum over support, ignoring NAs (C is #valid cells in support)
      out[] <- S
    } else {
      # only valid if the ENTIRE **box** has no NAs
      full_ok <- (round(Cc) >= scope_mass)
      out[full_ok] <- S[full_ok]
    }
  } else { # mean
    if (na.rm) {
      # average over valid cells in support
      Cz <- Cc
      Cz[Cz < tol] <- NA_real_
      out[] <- S / Cz
    } else {
      # mean only if ENTIRE **box** has no NAs; divide by sum of (actual) kernel
      ksum <- sum(kernel)
      full_ok <- (round(Cc) >= scope_mass)
      out[full_ok] <- S[full_ok] / ksum
    }
  }
  
  # terra's na.policy="omit": keep NA centers as NA
  if (na.policy == "omit") out[!center_valid] <- NA_real_
  
  # Keep orientation to match how fastfocal() writes values back
  t(out)
}

#' Next 5-smooth (2·3·5) length for FFT padding
#'
#' Rounds \code{n} up to the next integer whose prime factors are only
#' 2, 3, or 5. This size is typically fast for common FFT implementations.
#'
#' @param n Integer length.
#' @param factors Integer vector of allowed prime factors (default 2, 3, 5).
#'
#' @return Integer >= \code{n} that is "5-smooth".
#' @keywords internal
next_fast_len <- function(n, factors = c(2L,3L,5L)) {
  if (n <= 2L) return(2L)
  k <- as.integer(n)
  is_smooth <- function(x) {
    y <- x
    for (p in factors) {
      while ((y %% p) == 0L) y <- y %/% p
    }
    y == 1L
  }
  while (!is_smooth(k)) k <- k + 1L
  k
}
