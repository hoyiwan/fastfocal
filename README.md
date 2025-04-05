
# fastfocal <img src="https://img.shields.io/badge/license-GPL--3-blue.svg" alt="License: GPL-3" align="right"/>

`fastfocal` is a lightweight, high-performance R package for extracting raster data across spatial scales. It provides both **moving window (focal)** statistics and **point-based or buffer-based extraction** from raster layers — optimized for multi-scale modeling, ecological analysis, and machine learning workflows.

---

## 🚀 Features

- 🔁 Fast **moving window (focal)** raster operations across one or more spatial scales
- 📍 **Point-based or buffer-based** raster value extraction
- 📊 Supports multiple statistics: `"mean"` (default), `"sum"`, `"min"`, `"max"`, `"sd"`, and more
- ⚡ Built for speed with C++ backends
- 🧱 Compatible with both `terra::SpatRaster` and vector inputs from `sf` or `SpatVector`

---

## 📦 Installation

Install the development version directly from GitHub:

```r
# install.packages("devtools")  # if not installed
devtools::install_github("hoyiwan/fastfocal")
```

---

## 🛠️ Example Usage

### 🔁 Focal (Moving Window) Statistics with `fastfocal()`

```r
library(terra)
library(fastfocal)

# Load or create a raster with a projected CRS (units in meters)
r <- rast(matrix(runif(100), nrow = 10), extent = ext(0, 1000, 0, 1000))
crs(r) <- "EPSG:32610"  # UTM projection (meters)

# 1. Apply focal mean using a circular window with 500 m radius
out_500 <- fastfocal(r, window_size = 500, stat = "mean", shape = "circle")

# 2. Apply focal standard deviation using circular windows with radii of 100, 500, and 1000 meters
out_multi <- fastfocal(r, window_size = c(100, 500, 1000), stat = "sd", shape = "circle")

# Plot one of the results
plot(out_multi[[2]])  # e.g., result for 500m
```

> 💡 `window_size` is in **map units** (e.g., meters). Provide multiple values to apply focal stats at multiple spatial scales. `shape = "circle"` by default, can also use `"square"`.

---

## 📍 Extraction with `fastextract()`

```r
# Load point data (SpatVector or sf)
pts <- vect("points.shp")

# 1. Extract raster values directly at each point (default)
out_points <- fastextract(raster_stack, pts, stat = "mean")

# 2. Extract values within buffer zones at specified scales
out_buffers <- fastextract(raster_stack, pts, scales = c(500, 1000), stat = "mean")
```

> 💡 By default, `fastextract()` extracts values at point locations. Set `scales` to extract using buffer zones (in map units).

---

## 📄 License

This package is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

---

## 🤝 Contributing

Bug reports, feature suggestions, and contributions are welcome!  
Please open an [Issue](https://github.com/hoyiwan/fastfocal/issues) or submit a pull request.

---

## 📫 Author

**Ho Yi Wan**  
Associate Professor, University of Florida  
📧 hoyiwan@gmail.com
