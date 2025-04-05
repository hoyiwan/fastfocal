
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

### 🔁 Focal (Moving Window) Statistics

```r
library(terra)
library(fastfocal)

# Load or create a raster
r <- rast(matrix(1:100, nrow = 10))

# Apply focal mean using a 3x3 window
out1 <- fast_focal(r, window_size = 3, stat = "mean")

# Run across multiple scales (e.g., 3x3, 5x5, 7x7)
scales <- c(3, 5, 7)
out_stack <- focal_multi(r, scales = scales, stat = "sd")

# Plot result
plot(out1)
```

---

### 📍 Extraction with `fastextract()`

```r
# Load point data (SpatVector or sf)
pts <- vect("points.shp")

# 1. Extract raster values directly at each point (default)
out_points <- fastextract(raster_stack, pts, stat = "mean")

# 2. Extract values within buffer zones at specified scales
out_buffers <- fastextract(raster_stack, pts, scales = c(500, 1000), stat = "mean")
```

> 💡 `scales` defines buffer radii in map units. If omitted, extraction is done at point locations only.

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
