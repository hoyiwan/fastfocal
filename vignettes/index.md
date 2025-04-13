---
title: "fastfocal: Fast Multi-scale Raster Extraction and Moving Window Analysis with Fast Fourier Transform (FFT)"
subtitle: "Version 0.1.0 (2025)"
author: "Ho Yi Wan"
output: rmarkdown::html_vignette
bibliography: references.bib
vignette: >
  %\VignetteIndexEntry{fastfocal: Fast Multi-scale Raster Extraction and Moving Window Analysis with Fast Fourier Transform (FFT)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## What is `fastfocal`?

Welcome! I created the `fastfocal` R package (Wan 2025) because I do a lot of multi-scale ecological modeling, which requires applying focal statistics at various spatial scales. As rasters and kernels grow larger, traditional methods quickly become painfully slow — and I couldn’t find a solution that was both fast and flexible.

So I built one.

If you’ve landed here, chances are you’ve run into the same problem: you need to perform spatial smoothing or buffer-based extraction across large rasters or at multiple scales, and existing tools just aren’t fast or friendly enough. That’s exactly the need `fastfocal` was built to solve.

The `fastfocal` R package is designed for **high-performance focal and buffer-based raster processing**. It provides intuitive and fast **multi-layer aware** moving window operations and point-based extraction, with support for a rich variety of **kernel shapes**. Crucially, it automatically switches to **Fast Fourier Transform (FFT)** to optimize performance, dramatically accelerating large-kernel operations.Built on the `terra` framework, `fastfocal` works seamlessly with `SpatRaster` and `SpatVector` objects. It streamlines buffer extraction into a single-step function with **on-the-fly kernel generation**, so users can focus on analysis rather than boilerplate setup.

This vignette walks through the core functionality through simulated ecological raster layers and step-by-step examples. We highlight how `fastfocal()` and `fastextract()` offer smooth multi-scale analysis with flexible kernel definitions and performance-tuned design.

---

## 🌟 Key Features

- **Multi-layer raster support**: Apply focal operations across stacked rasters without manual iteration.
- **Flexible kernels**: Choose from 10+ window types (e.g., circle, Gaussian, logistic, quartic).
- **Map unit radius (in meters)**: Specify smoothing or extraction radius directly in meters.
- **Point-based extraction**: Extract raster summaries at points with optional buffer windows.
- **FFT backend**: Automatically switches to fast Fourier transform for large kernels.
- **Identical geometry output**: Outputs maintain resolution, extent, and CRS of the input.

---

## Installation and load packages


``` r
# run:
# install.packages("devtools")
# devtools::install_github("hoyiwan/fastfocal")

# load packages
library(fastfocal)
library(terra)

```

## 🌿 Simulating Rasters for This Vignette

To get started, we'll generate three example rasters with values resembling common ecological variables. In practice, users can replace these with their own data — `fastfocal` works seamlessly with any `SpatRaster` layers.


``` r
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
```

<div class="figure">
<img src="figure/simulate-coarse-fine-1.png" alt="plot of chunk simulate-coarse-fine" width="100%" />
<p class="caption">plot of chunk simulate-coarse-fine</p>
</div>

## 📍 Point & Buffer-Based Extraction Example

Here we'll use `fastextract` to extract values at a few points, then show how easy it is to pull in stats from buffers of different sizes — all in one step.

For this vignette, we simulate 3 points as an example. In practice, users can simply provide their own points or polygons as a `SpatVector`.


``` r
# simulate points
pts <- vect(matrix(c(500,500, 1500,1500, 2500,2500), ncol=2, byrow=TRUE), type="points", crs=crs(r))

# plot
par(mfrow = c(2, 2))
plot(r[[1]], main = "NDVI"); plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[2]], main = "Slope"); plot(pts, add = TRUE, pch = 20, cex = 2)
plot(r[[3]], main = "Temp");  plot(pts, add = TRUE, pch = 20, cex = 2)
par(mfrow = c(1, 1))
```

<div class="figure">
<img src="figure/simulate-points-1.png" alt="plot of chunk simulate-points" width="100%" />
<p class="caption">plot of chunk simulate-points</p>
</div>


``` r
# extract
fastextract(r, pts)  # defaults to d = 0 → extract at the point
#>        ndvi     slope     temp
#> 1 0.6658126 19.670246 14.81806
#> 2 0.6949464  9.219483 18.15735
#> 3 0.6780907 14.060074 21.59218
fastextract(r, pts, d = c(0, 100, 500))  # multiple distances
#>   scale_m      ndvi     slope     temp
#> 1       0 0.6658126 19.670246 14.81806
#> 2       0 0.6949464  9.219483 18.15735
#> 3       0 0.6780907 14.060074 21.59218
#> 4     100 0.6360472 18.382909 14.71503
#> 5     100 0.6728916  9.907261 18.03447
#> 6     100 0.6570599 14.158440 21.57079
#> 7     500 0.4974165 10.071724 14.27471
#> 8     500 0.5927882 10.990221 18.42000
#> 9     500 0.4596431 12.622186 21.19114
fastextract(r, pts, d = c(0, 100, 500, 1000), w = "gaussian", fun = "sd")
#>    scale_m       ndvi     slope      temp
#> 1        0         NA        NA        NA
#> 2        0         NA        NA        NA
#> 3        0         NA        NA        NA
#> 4      100 0.07220881 2.9682259 0.3835902
#> 5      100 0.05762610 1.7198168 0.2811265
#> 6      100 0.05638283 0.6116993 0.1554291
#> 7      500 0.16700735 4.0240024 1.2851331
#> 8      500 0.13796024 4.0145802 1.0747937
#> 9      500 0.14607090 3.0540960 1.0370660
#> 10    1000 0.16743288 4.7700249 2.0488525
#> 11    1000 0.14806001 4.3093104 1.9286310
#> 12    1000 0.16439540 3.5861606 1.7249144
```

## 🔁 Moving-window analysis

Finally, the fun part! 

We'll apply `fastfocal()` to our rasters and run moving-window stats across the whole stack using different kernels and summary functions. You’ll see how easy it is to explore multi-layer patterns with just a couple lines of code.


``` r
# Compute focal means using a circular window (300m radius)
f_circ <- fastfocal(r, d = 300, w = "circle", fun = "mean")
#> Using FFT backend\nUsing FFT backend\nUsing FFT backend\n

# Compute focal means using a Gaussian window (300m radius)
f_gaus <- fastfocal(r, d = 300, w = "gaussian", fun = "mean")
#> Using FFT backend\nUsing FFT backend\nUsing FFT backend\n

# plot
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
plot(r[[1]], main = "NDVI");  plot(f_circ[[1]], main = "NDVI - Circle");  plot(f_gaus[[1]], main = "NDVI - Gaussian")
plot(r[[2]], main = "Slope"); plot(f_circ[[2]], main = "Slope - Circle"); plot(f_gaus[[2]], main = "Slope - Gaussian")
plot(r[[3]], main = "Temp");  plot(f_circ[[3]], main = "Temp - Circle");  plot(f_gaus[[3]], main = "Temp - Gaussian")
```

<div class="figure">
<img src="figure/focal-1.png" alt="plot of chunk focal" width="100%" />
<p class="caption">plot of chunk focal</p>
</div>

``` r
par(mfrow = c(1, 1))
```

Let's see how different stats look like on the same layer!


``` r
# Compute different focal stats while holding the scale constant
f_mean <- fastfocal(r, d = 300, fun = "mean", na.rm = TRUE)
#> Using FFT backend\nUsing FFT backend\nUsing FFT backend\n
f_sd <- fastfocal(r, d = 300, fun = "sd", na.rm = TRUE)
#> Using terra::focal backend\nUsing terra::focal backend\nUsing terra::focal backend\n
f_median <- fastfocal(r, d = 300, fun = "median", na.rm = TRUE)
#> Using terra::focal backend\nUsing terra::focal backend\nUsing terra::focal backend\n

# plot
par(mfrow = c(2, 2))
plot(r[[1]], main = "NDVI")
plot(f_mean[[1]], main = "NDVI - Mean (300m)")
plot(f_sd[[1]], main = "NDVI - SD (300m)")
plot(f_median[[1]], main = "NDVI - Median (300m)")
```

<div class="figure">
<img src="figure/multi-stat-1.png" alt="plot of chunk multi-stat" width="100%" />
<p class="caption">plot of chunk multi-stat</p>
</div>

``` r
par(mfrow = c(1, 1))
```

## 🪟 Visualize Supported Kernels / Window Types

Curious what each kernel actually looks like? This section visualizes all the built-in window types supported by `fastfocal_weights()`. From sharp-edged circles to smooth tapers and heavy-tailed decay functions, each one defines how values are weighted across space during focal operations. Seeing them side-by-side makes it easier to choose the right kernel for your analysis.

The plot below shows each kernel shape, centered on a dummy raster with 30m resolution. The values are padded for contrast and normalized for consistent color scaling.


``` r
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
```

<div class="figure">
<img src="figure/visualize-kernels-1.png" alt="plot of chunk visualize-kernels" width="100%" />
<p class="caption">plot of chunk visualize-kernels</p>
</div>

``` r
par(mfrow = c(1, 1))
```

## 🧠 Notes

- Multi-layer support is now fully built into `fastfocal()`
- Input can be a single-layer raster or a stack
- Output preserves layer names and returns a `SpatRaster` with the same geometry
- Kernel size `d` is in **map units** (e.g., meters)
- Kernels were visualized using a center-pixel matrix to show their spatial shape

## 📋 Session Info


``` r
sessionInfo()
#> R version 4.4.3 (2025-02-28 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26100)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/Los_Angeles
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods  
#> [7] base     
#> 
#> other attached packages:
#> [1] fastfocal_0.1.0 terra_1.8-42   
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.9        digest_0.6.37     magrittr_2.0.3   
#>  [4] evaluate_1.0.3    pkgload_1.4.0     fastmap_1.2.0    
#>  [7] jsonlite_1.9.1    rprojroot_2.0.4   processx_3.8.6   
#> [10] pkgbuild_1.4.6    sessioninfo_1.2.3 rcmdcheck_1.4.0  
#> [13] urlchecker_1.0.1  ps_1.9.0          promises_1.3.2   
#> [16] purrr_1.0.4       codetools_0.2-20  jquerylib_0.1.4  
#> [19] cli_3.6.4         shiny_1.10.0      rlang_1.1.5      
#> [22] ellipsis_0.3.2    remotes_2.5.0     withr_3.0.2      
#> [25] cachem_1.1.0      yaml_2.3.10       devtools_2.4.5   
#> [28] tools_4.4.3       memoise_2.0.1     httpuv_1.6.15    
#> [31] curl_6.2.1        vctrs_0.6.5       R6_2.6.1         
#> [34] mime_0.12         lifecycle_1.0.4   fs_1.6.5         
#> [37] htmlwidgets_1.6.4 usethis_3.1.0     xopen_1.0.1      
#> [40] miniUI_0.1.1.1    pkgconfig_2.0.3   desc_1.4.3       
#> [43] callr_3.7.6       pillar_1.10.1     bslib_0.9.0      
#> [46] later_1.4.1       glue_1.8.0        profvis_0.4.0    
#> [49] Rcpp_1.0.14       xfun_0.51         tibble_3.2.1     
#> [52] rstudioapi_0.17.1 knitr_1.49        xtable_1.8-4     
#> [55] htmltools_0.5.8.1 rmarkdown_2.29    compiler_4.4.3   
#> [58] prettyunits_1.2.0
```

## 📚 Citation

To cite the `fastfocal` package in publications, use:

Wan, H.Y. (2025). fastfocal: Fast Multi-scale Raster Extraction and Moving Window Analysis with FFT. R package version X.X.X. https://github.com/hoyiwan/fastfocal


``` r
citation("fastfocal")
#> To cite package 'fastfocal' in publications use:
#> 
#>   Wan H (2025). _fastfocal: Fast Multi-Scale Focal and
#>   Extraction Functions for Raster Data_. R package version
#>   0.1.0, <https://github.com/hoyiwan/fastfocal>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {fastfocal: Fast Multi-Scale Focal and Extraction Functions for Raster Data},
#>     author = {Ho Yi Wan},
#>     year = {2025},
#>     note = {R package version 0.1.0},
#>     url = {https://github.com/hoyiwan/fastfocal},
#>   }
```
