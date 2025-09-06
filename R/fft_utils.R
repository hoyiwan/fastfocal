#' Next 5-smooth length (internal)
#'
#' Returns the smallest n' >= n whose prime factors are only 2, 3, and 5.
#' Useful to stabilize FFT performance by avoiding large-prime sizes.
#'
#' @param n integer length.
#' @keywords internal
#' @noRd
next_fast_len <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) return(1L)
  is_fast <- function(x) {
    y <- x
    for (p in c(2L, 3L, 5L)) while (y %% p == 0L) y <- y %/% p
    y == 1L
  }
  res <- n
  while (!is_fast(res)) res <- res + 1L
  res
}

# Internal: 2D FFT convolution with zero padding and valid-size crop
# Ensures a real-valued (numeric) result via Re(iFFT).
# @keywords internal
# @noRd
conv2_fft_ <- function(x, k, pad_nr, pad_nc) {
  nr <- nrow(x); nc <- ncol(x)
  kr <- nrow(k); kc <- ncol(k)
  want_nr <- nr + kr - 1L
  want_nc <- nc + kc - 1L
  
  X <- matrix(0, nrow = pad_nr, ncol = pad_nc)
  K <- matrix(0, nrow = pad_nr, ncol = pad_nc)
  X[1:nr, 1:nc] <- x
  K[1:kr, 1:kc] <- k
  
  FX <- stats::fft(X)
  FK <- stats::fft(K)
  Cc <- stats::fft(FX * FK, inverse = TRUE) / (pad_nr * pad_nc)
  
  # take the real part to avoid complex tiny imaginaries
  C <- Re(Cc)
  
  ro <- floor(kr / 2); co <- floor(kc / 2)
  r1 <- ro + 1L; r2 <- ro + nr
  c1 <- co + 1L; c2 <- co + nc
  C[r1:r2, c1:c2, drop = FALSE]
}

#' Masked FFT convolution with NA handling (internal)
#'
#' NA semantics:
#' - na.rm = TRUE: weighted mean/sum over available cells (kernel support only).
#' - na.rm = FALSE, na.policy = "omit": NA if any NA inside the full kr x kc box.
#' - na.rm = FALSE, na.policy = "all" : NA only if the entire box is missing.
#'
#' Returns a matrix aligned to the input and transposed at the end to match
#' terra::values orientation used in fastfocal().
#'
#' @param x numeric matrix (rows x cols).
#' @param kernel numeric matrix (unnormalized weights).
#' @param fun character, "mean" or "sum".
#' @param na.rm logical.
#' @param na.policy character, "omit" or "all".
#' @param pad character, "none" or "auto".
#' @keywords internal
#' @noRd
#' @importFrom stats fft
fft_convolve_masked <- function(x, kernel,
                                fun = "mean",
                                na.rm = TRUE,
                                na.policy = c("omit", "all"),
                                pad = c("none", "auto")) {
  stopifnot(is.matrix(x), is.matrix(kernel))
  na.policy <- match.arg(na.policy)
  pad <- match.arg(pad)
  
  nr <- nrow(x); nc <- ncol(x)
  kr <- nrow(kernel); kc <- ncol(kernel)
  
  # Valid mask and zero-filled values
  mask <- is.finite(x)
  xz <- x
  xz[!mask] <- 0
  
  # Padded sizes (for linear convolution)
  want_nr <- nr + kr - 1L
  want_nc <- nc + kc - 1L
  pad_nr  <- if (pad == "auto") next_fast_len(want_nr) else want_nr
  pad_nc  <- if (pad == "auto") next_fast_len(want_nc) else want_nc
  
  # Convolutions:
  #  - conv_val: values with the actual kernel
  #  - conv_w  : sum of kernel weights at valid cells (for mean when na.rm=TRUE)
  #  - conv_box: count of valid cells in a full kr x kc box (for gating when na.rm=FALSE)
  conv_val <- conv2_fft_(xz, kernel, pad_nr, pad_nc)
  
  support_mask <- (kernel != 0) * 1
  conv_w <- conv2_fft_(mask * 1, support_mask, pad_nr, pad_nc)
  
  box <- matrix(1, nrow = kr, ncol = kc)
  conv_box <- conv2_fft_(mask * 1, box, pad_nr, pad_nc)
  full_box <- kr * kc
  
  out <- conv_val
  
  if (identical(fun, "mean")) {
    # na.rm = TRUE: divide by weight sum on support
    nz <- conv_w != 0
    out[!nz] <- NA_real_
    out[nz] <- out[nz] / conv_w[nz]
  } else if (!identical(fun, "sum")) {
    stop("fft_convolve_masked supports only 'mean' or 'sum'.")
  }
  
  # Center NA handling
  if (na.policy == "omit") {
    out[!mask] <- NA_real_
  }
  
  # Gating for na.rm = FALSE
  if (!na.rm) {
    if (na.policy == "omit") {
      # Require the full box to be valid
      out[conv_box < full_box] <- NA_real_
      if (identical(fun, "mean")) {
        # When the box is fully valid, divide by sum(kernel) to match terra
        ksum <- sum(kernel)
        ok <- conv_box >= full_box
        out[ok] <- out[ok] / ksum
      }
    } else { # na.policy == "all"
      # Only require at least one valid cell in the box
      out[conv_box == 0] <- NA_real_
      if (identical(fun, "mean")) {
        # Divide by sum(kernel) for mean when not removing NAs
        ksum <- sum(kernel)
        nz <- conv_box != 0
        out[nz] <- out[nz] / ksum
      }
    }
  }
  
  # Transpose to match terra::values orientation
  t(out)
}
