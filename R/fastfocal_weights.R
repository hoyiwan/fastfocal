#' Generate weight matrix for focal operations using map units
#'
#' Builds an UNnormalized kernel from map units. Circle uses a center-distance
#' rule (include a cell if its center is within radius `d`); other kernels follow
#' the same definitions as before. For exact equality with terra's discrete
#' circle, pass `w = terra::focalMat(x, d, type = "circle")` into `fastfocal()`.
#'
#' @param x SpatRaster (for resolution).
#' @param d numeric. Radius in map units.
#' @param w character. "rectangle","circle","circular","gaussian","Gauss",
#'   "pareto","idw","exponential","triangular","cosine","logistic","cauchy",
#'   "quartic","epanechnikov".
#' @param normalize logical. Normalize weights to sum to 1 (default TRUE).
#' @param plot logical. If TRUE, plots the kernel.
#' @return numeric matrix of weights.
#' @export
#' @importFrom terra res
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors
fastfocal_weights <- function(x, d, w = "circle", normalize = TRUE, plot = FALSE) {
  if (!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
  stopifnot(d > 0)
  
  res_val <- terra::res(x)[1]  # assumes square pixels
  cell_radius <- ceiling(d / res_val)
  size <- 2L * cell_radius + 1L
  center <- ceiling(size / 2)
  
  mat <- matrix(0, nrow = size, ncol = size)
  
  # Precompute center distances for speed
  ii <- matrix(rep(seq_len(size), size), nrow = size)
  jj <- t(ii)
  dx <- (ii - center) * res_val
  dy <- (jj - center) * res_val
  dist <- sqrt(dx*dx + dy*dy)
  
  if (w %in% c("rectangle")) {
    mat[,] <- 1
  } else if (w %in% c("circle","circular")) {
    mat[dist <= d] <- 1
  } else if (w == "gaussian" || w == "Gauss") {
    sigma <- d / 3
    mat <- exp(-(dist^2) / (2 * sigma^2))
  } else if (w == "pareto") {
    mat[dist <= d] <- (1 + dist[dist <= d] / (d / 3))^-2
  } else if (w == "idw") {
    mat[dist == 0] <- 1
    inside <- dist > 0 & dist <= d
    mat[inside] <- 1 / (dist[inside] / (d / 3) + 1)^2
  } else if (w == "exponential") {
    inside <- dist <= d
    mat[inside] <- exp(-dist[inside] / (d / 3))
  } else if (w == "triangular") {
    inside <- dist <= d
    mat[inside] <- 1 - dist[inside] / d
  } else if (w == "cosine") {
    inside <- dist <= d
    mat[inside] <- cos(pi * dist[inside] / (2 * d))^2
  } else if (w == "logistic") {
    steepness <- 10 / d
    mat <- 1 / (1 + exp(steepness * (dist - d / 2)))
  } else if (w == "cauchy") {
    inside <- dist <= d
    mat[inside] <- 1 / (1 + (dist[inside] / (d / 3))^2)
  } else if (w == "quartic") {
    inside <- dist <= d
    mat[inside] <- (1 - (dist[inside] / d)^2)^2
  } else if (w == "epanechnikov") {
    inside <- dist <= d
    mat[inside] <- 1 - (dist[inside] / d)^2
  } else {
    stop("Unsupported kernel type: ", w)
  }
  
  if (normalize) {
    s <- sum(mat)
    if (s > 0) mat <- mat / s
  }
  
  if (plot) {
    graphics::image(mat, main = paste("fastfocal_weights:", w),
                    col = grDevices::hcl.colors(20, "YlOrRd", rev = TRUE))
  }
  
  mat
}
