## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  fig.width = 6,
  fig.height = 4,
  out.width = "100%"
)
library(terra)
library(fastfocal)

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("devtools")
# # devtools::install_github("hoyiwan/fastfocal")

## -----------------------------------------------------------------------------
library(fastfocal)
library(terra)

## ----simulate-coarse-fine-----------------------------------------------------
set.seed(888)
ndvi_coarse  <- rast(matrix(runif(100, 0.1, 0.9), 10, 10))            # NDVI (0.1 - 0.9)
slope_coarse <- rast(matrix(rbeta(100, 2, 5) * 40, 10, 10))            # slope (degrees)
temp_coarse  <- rast(matrix(seq(12, 24, length.out = 100) + rnorm(100, 0, 1), 10, 10))  # temperature with gradient

# Disaggregate to finer resolution
ndvi  <- disagg(ndvi_coarse, fact = 10, method = "bilinear")
slope <- disagg(slope_coarse, fact = 10, method = "bilinear")
temp  <- disagg(temp_coarse, fact = 10, method = "bilinear")

# Set a 3 km by 3 km extent for all layers (projected units, e.g., meters)
ext(ndvi) <- ext(0, 3000, 0, 3000)
ext(slope) <- ext(ndvi)
ext(temp) <- ext(ndvi)
crs(ndvi) <- "EPSG:32611"  # UTM Zone 11N
crs(slope) <- crs(ndvi)
crs(temp) <- crs(ndvi)

# Combine into a multi-layer SpatRaster
r <- c(ndvi, slope, temp)
names(r) <- c("ndvi", "slope", "temp")

plot(r)

## ----simulate-points----------------------------------------------------------
# Simulate some points
pts <- vect(matrix(c(500,500, 1500,1500, 2500,2500), ncol=2, byrow=TRUE), type="points", crs=crs(r))

par(mfrow = c(2, 2))
plot(r[[1]], main = "NDVI");  plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[2]], main = "Slope"); plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[3]], main = "Temp");  plot(pts, add = TRUE, pch = 20, cex = 2)
par(mfrow = c(1, 1))

## ----extract------------------------------------------------------------------
fastextract(r, pts)                               # default d = 0 (extract at the point)
fastextract(r, pts, d = c(0, 100, 500))           # multiple distances
fastextract(r, pts, d = c(0, 100, 500, 1000), w = "gaussian", fun = "sd")

## ----focal--------------------------------------------------------------------
# Focal means using a circular window (300 m radius)
f_circ <- fastfocal(r, d = 300, w = "circle",   fun = "mean")

# Focal means using a gaussian window (300 m radius)
f_gaus <- fastfocal(r, d = 300, w = "gaussian", fun = "mean")

par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
plot(r[[1]], main = "NDVI");  plot(f_circ[[1]], main = "NDVI - Circle");  plot(f_gaus[[1]], main = "NDVI - Gaussian")
plot(r[[2]], main = "Slope"); plot(f_circ[[2]], main = "Slope - Circle"); plot(f_gaus[[2]], main = "Slope - Gaussian")
plot(r[[3]], main = "Temp");  plot(f_circ[[3]], main = "Temp - Circle");  plot(f_gaus[[3]], main = "Temp - Gaussian")
par(mfrow = c(1, 1))

## ----multi-stat---------------------------------------------------------------
f_mean   <- fastfocal(r, d = 300, fun = "mean",   na.rm = TRUE)
f_sd     <- fastfocal(r, d = 300, fun = "sd",     na.rm = TRUE)
f_median <- fastfocal(r, d = 300, fun = "median", na.rm = TRUE)

par(mfrow = c(2, 2))
plot(r[[1]],       main = "NDVI")
plot(f_mean[[1]],  main = "NDVI - Mean (300 m)")
plot(f_sd[[1]],    main = "NDVI - SD (300 m)")
plot(f_median[[1]],main = "NDVI - Median (300 m)")
par(mfrow = c(1, 1))

## ----visualize-kernels--------------------------------------------------------
center_r <- rast(ext(0, 90, 0, 90), resolution = 30, crs = "EPSG:32611")
kernel_types <- c("circle", "rectangle", "gaussian", "pareto", "idw", "exponential",
                  "triangular", "cosine", "logistic", "cauchy", "quartic", "epanechnikov")

pad_kernel <- function(k, pad = 2) {
  nr <- nrow(k); nc <- ncol(k)
  out <- matrix(0, nrow = nr + 2 * pad, ncol = nc + 2 * pad)
  out[(pad + 1):(pad + nr), (pad + 1):(pad + nc)] <- k
  out
}

par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))
for (w in kernel_types) {
  k_raw <- fastfocal_weights(center_r, d = 150, w = w, normalize = FALSE)  # ~11x11 at 30 m
  k_pad <- pad_kernel(k_raw, pad = 2)
  k <- k_pad / max(k_pad, na.rm = TRUE)
  image(k, main = w, col = hcl.colors(20, "YlOrRd", rev = TRUE), zlim = c(0, 1), asp = 1, cex.main = 1.2)
}
par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
citation("fastfocal")

