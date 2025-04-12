## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE,
  fig.width = 6,
  fig.height = 4,
  out.width="100%"
)
library(terra)
library(fastfocal)

## -----------------------------------------------------------------------------
# run:
# install.packages("devtools")
# devtools::install_github("hoyiwan/fastfocal")

# load packages
library(fastfocal)
library(terra)


## ----simulate-coarse-fine-----------------------------------------------------
# create random rasters
set.seed(888)
ndvi_coarse  <- rast(matrix(runif(100, 0.1, 0.9), 10, 10)) # Simulate NDVI values (0.1 - 0.9)
slope_coarse <- rast(matrix(rbeta(100, 2, 5) * 40, 10, 10)) # simulate slope values (0 - 40 degrees)
temp_coarse  <- rast(matrix(seq(12, 24, length.out = 100) + rnorm(100, 0, 1), 10, 10)) # simulate temperature values with a gentle gradient + random noise

# Disaggregate each raster to make them look more natural
ndvi  <- disagg(ndvi_coarse, fact = 10, method = "bilinear")
slope <- disagg(slope_coarse, fact = 10, method = "bilinear")
temp  <- disagg(temp_coarse, fact = 10, method = "bilinear")

# Set extent to 3km × 3km for all layers (in projected units, e.g., meters)
ext(ndvi) <- ext(0, 3000, 0, 3000)
ext(slope) <- ext(ndvi)
ext(temp) <- ext(ndvi)
crs(ndvi) <- "EPSG:32611" # Assign UTM coordinate system (Zone 11N) to all layers
crs(slope) <- crs(ndvi)
crs(temp) <- crs(ndvi)

# Combine layers into a multi-layer SpatRaster stack and name them
r <- c(ndvi, slope, temp)
names(r) <- c("ndvi", "slope", "temp")

# plot
plot(r)

## ----simulate-points----------------------------------------------------------
# simulate points
pts <- vect(matrix(c(500,500, 1500,1500, 2500,2500), ncol=2, byrow=TRUE), type="points", crs=crs(r))

# plot
par(mfrow = c(2, 2))
plot(r[[1]], main = "NDVI"); plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[2]], main = "Slope"); plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[3]], main = "Temp");  plot(pts, add = TRUE, pch = 20, cex = 2)
par(mfrow = c(1, 1))

## ----extract------------------------------------------------------------------
# extract
fastextract(r, pts)  # defaults to d = 0 → extract at the point
fastextract(r, pts, d = c(0, 100, 500))  # multiple distances
fastextract(r, pts, d = c(0, 100, 500, 1000), w = "gaussian", fun = "sd")

## ----focal--------------------------------------------------------------------
# Compute focal means using a circular window (300m radius)
f_circ <- fastfocal(r, d = 300, w = "circle", fun = "mean")

# Compute focal means using a Gaussian window (300m radius)
f_gaus <- fastfocal(r, d = 300, w = "gaussian", fun = "mean")

# plot
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
plot(r[[1]], main = "NDVI");  plot(f_circ[[1]], main = "NDVI - Circle");  plot(f_gaus[[1]], main = "NDVI - Gaussian")
plot(r[[2]], main = "Slope"); plot(f_circ[[2]], main = "Slope - Circle"); plot(f_gaus[[2]], main = "Slope - Gaussian")
plot(r[[3]], main = "Temp");  plot(f_circ[[3]], main = "Temp - Circle");  plot(f_gaus[[3]], main = "Temp - Gaussian")
par(mfrow = c(1, 1))

## ----multi-stat---------------------------------------------------------------
# Compute different focal stats while holding the scale constant
f_mean <- fastfocal(r, d = 300, fun = "mean", na.rm = TRUE)
f_sd <- fastfocal(r, d = 300, fun = "sd", na.rm = TRUE)
f_median <- fastfocal(r, d = 300, fun = "median", na.rm = TRUE)

# plot
par(mfrow = c(2, 2))
plot(r[[1]], main = "NDVI")
plot(f_mean[[1]], main = "NDVI - Mean (300m)")
plot(f_sd[[1]], main = "NDVI - SD (300m)")
plot(f_median[[1]], main = "NDVI - Median (300m)")
par(mfrow = c(1, 1))

## ----visualize-kernels--------------------------------------------------------
# simulate kernels
center_r <- rast(ext(0, 90, 0, 90), resolution = 30, crs = "EPSG:32611")
kernel_types <- c("circle", "rectangle", "gaussian", "pareto", "idw", "exponential",
                  "triangular", "cosine", "logistic", "cauchy", "quartic", "epanechnikov")

pad_kernel <- function(k, pad = 20) {
  nr <- nrow(k); nc <- ncol(k)
  out <- matrix(0, nrow = nr + 2 * pad, ncol = nc + 2 * pad)
  out[(pad + 1):(pad + nr), (pad + 1):(pad + nc)] <- k
  return(out)
}

# plot
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))
for (w in kernel_types) {
  k_raw <- fastfocal_weights(center_r, d = 1000, w = w, normalize = FALSE)
  k_pad <- pad_kernel(k_raw, pad = 20)
  k <- k_pad / max(k_pad, na.rm = TRUE)
  image(k, main = w, col = hcl.colors(20, "YlOrRd", rev = TRUE), zlim = c(0, 1), asp = 1, cex.main = 1.5)
}
par(mfrow = c(1, 1))

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
citation("fastfocal")


