#' Generate weight matrix for focal operations using map units
#'
#' Builds an unnormalized (or normalized) kernel from map units.
#' Circle uses a center-distance rule (include if center <= d).
#' **Gaussian interprets d as sigma in map units and truncates at 3 sigma**,
#' matching terra::focalMat(..., type = "Gauss").
#'
#' @param x SpatRaster (used for resolution; assumes square pixels).
#' @param d numeric. Radius in map units for most kernels; **sigma** in map units for "gaussian"/"Gauss".
#' @param w character. One of:
#'   "rectangle","circle","circular","gaussian","Gauss","pareto","idw",
#'   "exponential","triangular","cosine","logistic","cauchy","quartic","epanechnikov".
#' @param normalize logical. If TRUE (default), scale weights to sum to 1.
#' @param plot logical. If TRUE, plots the kernel.
#' @return numeric matrix of weights.
#' @export
#' 
#' @examples
#' # Small raster (resolution = 1 map unit)
#' r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
#'
#' # Circle: d is a radius in map units -> here cell_radius = 2 -> 5x5 kernel
#' Kc <- fastfocal_weights(r, d = 2, w = "circle", normalize = TRUE)
#' dim(Kc)            # 5 x 5
#' round(sum(Kc), 6)  # ~1
#'
#' # Gaussian: d is sigma in map units, truncated at 3 sigmas
#' Kg <- fastfocal_weights(r, d = 1, w = "gaussian", normalize = TRUE)
#' dim(Kg)            # 7 x 7 (since 2*ceil(3*sigma) + 1)
#' round(sum(Kg), 6)  # ~1
#'
#' # \donttest{
#' # Quick visualization (kept out of CRAN's main run)
#' fastfocal_weights(r, d = 2, w = "circle", normalize = TRUE, plot = TRUE)
#' # }

fastfocal_weights <- function(x, d, w = "circle", normalize = TRUE, plot = FALSE) {
  if (!inherits(x, "SpatRaster")) stop("x must be a SpatRaster")
  stopifnot(is.numeric(d), d > 0)
  
  res_val <- terra::res(x)[1]  # assumes square pixels
  
  # ---- window size (cells) ----
  if (w %in% c("gaussian", "Gauss")) {
    # terra-compatible: d is sigma; truncate at ~3*sigma
    sigma <- d
    cell_radius <- ceiling(3 * sigma / res_val)
  } else {
    # others (circle, etc): d is a radius
    cell_radius <- ceiling(d / res_val)
  }
  size   <- 2L * cell_radius + 1L
  center <- cell_radius + 1L
  
  # grid in map units
  ii <- matrix(rep(seq_len(size), size), nrow = size)
  jj <- t(ii)
  dx <- (ii - center) * res_val
  dy <- (jj - center) * res_val
  dist <- sqrt(dx * dx + dy * dy)
  
  mat <- matrix(0, nrow = size, ncol = size)
  
  if (w %in% c("rectangle")) {
    mat[] <- 1
    
  } else if (w %in% c("circle","circular")) {
    mat[dist <= d] <- 1
    
  } else if (w %in% c("gaussian","Gauss")) {
    # d is sigma; no hard cutoff other than the 3 sigma window we built
    sigma <- d
    mat <- exp(-(dist^2) / (2 * sigma^2))
    
  } else if (w == "pareto") {
    inside <- dist <= d
    mat[inside] <- (1 + dist[inside] / (d / 3))^-2
    
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
  
  if (isTRUE(plot)) {
    graphics::image(
      mat, main = paste("fastfocal_weights:", w),
      col = grDevices::hcl.colors(20, "YlOrRd", rev = TRUE)
    )
  }
  mat
}
