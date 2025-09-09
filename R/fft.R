# ---- fft.R ---------------------------------------------------------------
# Internal FFT helpers for fastfocal
# (not exported)

#' Next 5-smooth length (internal)
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

#' 2D FFT convolution with zero padding + valid-size crop (internal)
#' Aligns raw convolution back to input extents.
#' @keywords internal
#' @noRd
conv2_fft_ <- function(x, k, pad_nr, pad_nc) {
  nr <- nrow(x); nc <- ncol(x)
  kr <- nrow(k); kc <- ncol(k)
  
  X <- matrix(0, nrow = pad_nr, ncol = pad_nc)
  K <- matrix(0, nrow = pad_nr, ncol = pad_nc)
  X[1:nr, 1:nc] <- x
  K[1:kr, 1:kc] <- k
  
  FX <- stats::fft(X)
  FK <- stats::fft(K)
  Cc <- stats::fft(FX * FK, inverse = TRUE) / (pad_nr * pad_nc)
  
  C  <- Re(Cc)
  ro <- floor(kr / 2); co <- floor(kc / 2)
  r1 <- ro + 1L; r2 <- ro + nr
  c1 <- co + 1L; c2 <- co + nc
  C[r1:r2, c1:c2, drop = FALSE]
}

#' Masked FFT convolution with terra-parity NA semantics
#'
#' `x`, `kernel`: matrices in terra's `as.matrix(..., wide=TRUE)` orientation.
#' Returns **t(out)** (single transpose) because callers map via `values<-`.
#'
#' NA rules:
#' - `na.rm = TRUE`:
#'     * mean: weighted avg over available cells; NA iff no valid cell.
#'     * sum : sum over available cells; NA iff no valid cell.
#' - `na.rm = FALSE`:
#'     * omit: require **all overlapped** cells to be valid (edge-aware);
#'             mean divides by overlapped weight sum.
#'     * all : NA only if overlapped support has no valid cells;
#'             mean divides by overlapped weight sum.
#'
#' A small epsilon + integer rounding avoids FFT round-off speckles.
#'
#' @keywords internal
#' @noRd
fft_convolve_masked <- function(x, kernel,
                                fun = c("mean","sum"),
                                na.rm = TRUE,
                                na.policy = c("omit","all"),
                                pad = c("auto","none")) {
  stopifnot(is.matrix(x), is.matrix(kernel))
  fun       <- match.arg(fun)
  na.policy <- match.arg(na.policy)
  pad       <- match.arg(pad)
  
  nr <- nrow(x); nc <- ncol(x)
  kr <- nrow(kernel); kc <- ncol(kernel)
  full_size <- kr * kc
  ksum <- sum(kernel)
  
  # mask & padded sizes
  mask <- is.finite(x)
  xz   <- x; xz[!mask] <- 0
  
  want_nr <- nr + kr - 1L
  want_nc <- nc + kc - 1L
  pad_nr  <- if (pad == "auto") next_fast_len(want_nr) else want_nr
  pad_nc  <- if (pad == "auto") next_fast_len(want_nc) else want_nc
  
  # core convs (cropped to nr x nc)
  conv_val <- conv2_fft_(xz,        kernel,             pad_nr, pad_nc)  # sum(x*w) over valid cells
  conv_w   <- conv2_fft_(mask * 1,  kernel,             pad_nr, pad_nc)  # sum of weights over valid cells
  conv_box <- conv2_fft_(mask * 1,  matrix(1, kr, kc),  pad_nr, pad_nc)  # count of valid cells in overlap
  conv_sup <- conv2_fft_(matrix(1, nr, nc), matrix(1, kr, kc), pad_nr, pad_nc) # size of overlapped support
  
  # jitter guards
  eps        <- 1e-12
  conv_box_i <- round(conv_box)
  conv_sup_i <- round(conv_sup)
  has_w      <- (abs(conv_w) > eps)
  
  out <- conv_val
  
  if (isTRUE(na.rm)) {
    if (fun == "mean") {
      out[!has_w] <- NA_real_
      ok <- has_w
      out[ok] <- out[ok] / conv_w[ok]
    } else { # sum
      out[conv_box_i == 0L] <- NA_real_
    }
    
  } else { # na.rm == FALSE
    if (na.policy == "omit") {
      # terra rule: require a FULL window (not just the overlapped part) and all valid
      ok_full <- (conv_box_i == full_size)  # implies edge windows (conv_sup < full_size) are excluded
      out[!ok_full] <- NA_real_
      if (fun == "mean") {
        # full window => divide by constant kernel sum (terra parity)
        out[ok_full] <- out[ok_full] / ksum
      }
    } else { # na.policy == "all"
      # only all-NA overlap -> NA; mean divides by overlapped weight sum
      out[conv_box_i == 0L] <- NA_real_
      if (fun == "mean") {
        ok <- has_w
        out[ok] <- out[ok] / conv_w[ok]
      }
    }
  }
  
  # return in terra orientation (one transpose)
  t(out)
}
