#' Generate weight matrix for focal operations using map units
#'
#' @param x SpatRaster. Used to determine resolution.
#' @param d Radius in map units (e.g. meters).
#' @param w Type of weight. One of:
#'   "rectangle", "circle", "gaussian", "pareto", "idw", "exponential",
#'   "triangular", "cosine", "logistic", "cauchy", "quartic", "epanechnikov".
#' @param normalize Logical. Whether to normalize weights so they sum to 1. Default = TRUE.
#' @param plot Logical. If TRUE, plots the resulting kernel. Default = FALSE.
#'
#' @return A numeric matrix of focal weights.
#' @export
#'
#' @importFrom terra rast
#' 
#' @examples
#' library(terra)
#' x <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000)
#' values(x) <- runif(ncell(x))
#' fastfocal_weights(x, d = 300, w = "gaussian", plot = TRUE)
fastfocal_weights <- function(x, d, w = "circle", normalize = TRUE, plot = FALSE) {
  if (!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
  stopifnot(d > 0)
  
  res_val <- res(x)[1]  # assumes square pixels
  cell_radius <- ceiling(d / res_val)
  size <- 2 * cell_radius + 1
  center <- ceiling(size / 2)
  mat <- matrix(0, nrow = size, ncol = size)
  
  for (i in 1:size) {
    for (j in 1:size) {
      dx <- (i - center) * res_val
      dy <- (j - center) * res_val
      dist <- sqrt(dx^2 + dy^2)
      
      mat[i, j] <- switch(w,
                          "rectangle" = 1,
                          "circle" = ifelse(dist <= d, 1, 0),
                          "circular" = ifelse(dist <= d, 1, 0),
                          "gaussian" = {
                            sigma <- d / 3
                            exp(-(dist^2) / (2 * sigma^2))
                          },
                          "Gauss" = {
                            sigma <- d / 3
                            exp(-(dist^2) / (2 * sigma^2))
                          },
                          "pareto" = ifelse(dist <= d, (1 + dist / (d / 3))^-2, 0),
                          "idw" = ifelse(dist == 0, 1, ifelse(dist <= d, 1 / (dist / (d / 3) + 1)^2, 0)),
                          "exponential" = ifelse(dist <= d, exp(-dist / (d / 3)), 0),
                          "triangular" = ifelse(dist <= d, 1 - dist / d, 0),
                          "cosine" = ifelse(dist <= d, cos(pi * dist / (2 * d))^2, 0),
                          "logistic" = {
                            steepness <- 10 / d
                            1 / (1 + exp(steepness * (dist - d / 2)))
                          },
                          "cauchy" = ifelse(dist <= d, 1 / (1 + (dist / (d / 3))^2), 0),
                          "quartic" = ifelse(dist <= d, (1 - (dist / d)^2)^2, 0),
                          "epanechnikov" = ifelse(dist <= d, 1 - (dist / d)^2, 0),
                          stop("Unsupported kernel type")
      )
    }
  }
  
  if (normalize) {
    mat <- mat / sum(mat)
  }
  
  if (plot) {
    image(mat, main = paste("fastfocal_weights:", w), col = hcl.colors(20, "YlOrRd", rev = TRUE))
  }
  
  return(mat)
}