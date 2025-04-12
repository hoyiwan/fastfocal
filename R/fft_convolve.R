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
  x_mat <- as.matrix(x)
  
  # Replace NAs with 0 for FFT
  x_na_mask <- is.na(x_mat)
  x_mat[x_na_mask] <- 0
  
  # Zero-pad kernel to match x_mat dimensions
  nr <- nrow(x_mat) + nrow(kernel) - 1
  nc <- ncol(x_mat) + ncol(kernel) - 1
  
  pad_x <- matrix(0, nr, nc)
  pad_k <- matrix(0, nr, nc)
  pad_x[1:nrow(x_mat), 1:ncol(x_mat)] <- x_mat
  pad_k[1:nrow(kernel), 1:ncol(kernel)] <- kernel
  
  # FFT convolution
  fft_x <- fft(pad_x)
  fft_k <- fft(pad_k)
  conv <- Re(fft(fft_x * fft_k, inverse = TRUE)) / length(fft_x)
  
  # NA correction
  if (na.rm) {
    mask <- matrix(1, nrow = nrow(x_mat), ncol = ncol(x_mat))
    mask[x_na_mask] <- 0
    
    pad_m <- matrix(0, nr, nc)
    pad_m[1:nrow(mask), 1:ncol(mask)] <- mask
    norm <- Re(fft(fft(pad_m) * fft_k, inverse = TRUE)) / length(fft_x)
    
    norm[norm == 0] <- NA  # Avoid divide-by-zero
    if (fun == "mean") {
      conv <- conv / norm
    }
  } else {
    if (fun == "mean") {
      conv <- conv / sum(kernel)
    }
  }
  
  # Align and crop to match input dimensions
  row_off <- floor(nrow(kernel) / 2)
  col_off <- floor(ncol(kernel) / 2)
  result <- conv[
    (row_off + 1):(row_off + nrow(x_mat)),
    (col_off + 1):(col_off + ncol(x_mat))
  ]
  
  return(t(result))
}