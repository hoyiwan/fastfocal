# fastfocal

[![CRAN Status](https://www.r-pkg.org/badges/version/fastfocal)](https://CRAN.R-project.org/package=fastfocal)
[![R-CMD-check](https://github.com/hoyiwan/fastfocal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hoyiwan/fastfocal/actions)
[![DOI](https://zenodo.org/badge/961060307.svg)](https://doi.org/10.5281/zenodo.17074691)

**fastfocal: Fast Multi-scale Raster Extraction and Moving Window Analysis with Fast Fourier Transform (FFT) in R**

`fastfocal` provides high-performance, flexible raster smoothing and extraction functions in R using moving windows, buffer-based zones, and an auto-switching FFT backend for large kernels. It supports multiple focal statistics and allows users to work at multiple spatial scales with ease.

---

## Installation

Once accepted on CRAN, you will be able to install the stable release with:

```r
install.packages("fastfocal")
```

Until then, you can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hoyiwan/fastfocal")
```

---

## Overview

This package includes:

- `fastfocal()` — fast moving-window smoothing with support for mean, sum, min, max, sd, and median  
- `fastextract()` — fast extraction of raster values at point or buffer locations  
- `fastfocal_weights()` — utility for generating spatial weight matrices (circular, Gaussian, etc.)  
- Auto-switch backend to FFT for large windows to improve performance  
  - FFT backend currently supports **sum** and **mean**; other statistics use the `terra` backend  
- Native support for `terra::SpatRaster` and `terra::SpatVector` objects  

---

## Example Usage

```r
library(fastfocal)
library(terra)

# Create a dummy raster
r <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000)
values(r) <- runif(ncell(r))

# Apply fast focal smoothing with circular window of radius 300
smoothed <- fastfocal(r, d = 300, w = "circle", fun = "mean")

# Plot the result
plot(smoothed)
```

Weighted extraction at points:

```r
# Create SpatVector of points
pts <- vect(data.frame(x = c(500, 1500), y = c(500, 2500)), geom = c("x", "y"), crs = crs(r))

# Extract raster values in 500 m buffers around points
result <- fastextract(r, pts, d = 500, fun = "mean")
print(result)
```

---

## Vignettes

- [Introduction (Index)](https://hoyiwan.github.io/fastfocal/index.html)  
- [Benchmark: fastfocal vs terra::focal](https://hoyiwan.github.io/fastfocal/benchmark.html)  

You can also access them from R using:

```r
vignette("index", package = "fastfocal")
```

---

## License

This package is licensed under the **MIT License** (see the [LICENSE](LICENSE) file).

---

## Citation

If you use `fastfocal` in published work, please cite it as:

> Ho Yi Wan (2025). *fastfocal: A fast, energy-efficient R package for focal raster operations*. Version v0.1.1. Zenodo. https://doi.org/10.5281/zenodo.17074691

Or use the BibTeX entry:

```bibtex
@software{wan_fastfocal_2025,
  author       = {Ho Yi Wan},
  title        = {fastfocal: A fast, energy-efficient R package for focal raster operations},
  version      = {v0.1.3},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17074691},
  url          = {https://doi.org/10.5281/zenodo.17074691}
}
```

You can also run:

```r
citation("fastfocal")
```

---

## Author

**Ho Yi Wan**  
hoyiwan@gmail.com  

---

## Purpose

Built to support large-scale ecological analysis, high-performance raster processing, and reproducible landscape research workflows.
