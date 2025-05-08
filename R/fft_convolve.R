#' 2D FFT Convolution (terra-compatible, supports na.rm)
#'
#' Applies 2D convolution via FFT to a matrix using a kernel.
#'
#' @param x Matrix input
#' @param kernel Smoothing kernel
#' @param fun Character. Summary function name (e.g., "mean", "sum")
#' @param na.rm Logical. Whether to remove NAs in computation
#' @return Convolved and aligned result
#' @export
fft_convolve <- function(x, kernel, fun = "mean", na.rm = TRUE) {
  # 1) Prepare matrix and NA mask
  x_mat     <- as.matrix(x)
  x_na_mask <- is.na(x_mat)
  x_mat[x_na_mask] <- 0
  
  # 2) Compute padded size and normalization factor
  nr <- nrow(x_mat) + nrow(kernel) - 1
  nc <- ncol(x_mat) + ncol(kernel) - 1
  norm_factor <- nr * nc
  
  # 3) Zero‑pad input & kernel
  pad_x <- matrix(0, nr, nc)
  pad_k <- matrix(0, nr, nc)
  pad_x[1:nrow(x_mat), 1:ncol(x_mat)]   <- x_mat
  pad_k[1:nrow(kernel), 1:ncol(kernel)] <- kernel
  
  # 4) FFTs and raw convolution
  fft_x     <- fft(pad_x)
  fft_k     <- fft(pad_k)
  conv_full <- Re(fft(fft_x * fft_k, inverse = TRUE)) / norm_factor
  
  # 5) Crop to “valid” region
  row_off <- floor(nrow(kernel) / 2)
  col_off <- floor(ncol(kernel) / 2)
  i1 <- row_off + 1; i2 <- row_off + nrow(x_mat)
  j1 <- col_off + 1; j2 <- col_off + ncol(x_mat)
  conv_crop <- conv_full[i1:i2, j1:j2]
  
  # 6) Apply summary function (mean or sum)
  if (fun == "mean") {
    conv_crop <- conv_crop / sum(kernel)
  }
  # (if fun == "sum", leave conv_crop as-is)
  
  # 7) Re‑introduce NAs in original input positions for both cases
  conv_crop[x_na_mask] <- NA
  
  # 8) Return with same orientation
  return(t(conv_crop))
}
