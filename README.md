
# fastfocal

**fastfocal** is an R package for high-performance, multi-scale raster processing. It provides:

- ⚡️ Fast raster extraction at point locations or buffer radii
- 🔁 Fast moving window statistics (focal operations)
- 🧠 Support for multiple window shapes: `circular`, `rectangular`, `gaussian`
- 📐 Handles rasters with multiple layers
- 📊 Progress bar and flexible NA handling
- 💻 Powered by C++ for serious performance

---

## Installation

```r
# From GitHub
devtools::install_github("hoyiwan/fastfocal")
```

---

## Example: Point Extraction

```r
library(fastfocal)
library(terra)

# Create raster and points
r <- rast(matrix(1:100, 10, 10))
ext(r) <- c(0, 10, 0, 10)
pts <- vect(matrix(c(1,1, 5,5, 9,9), ncol=2, byrow=TRUE), type="points")

# Extract raster values at point and buffer scales
fastextract(r, pts, stat = "mean", scales = c(0, 2), na.rm = TRUE)
```

---

## Example: Moving Window Focal Statistics

```r
# Apply mean filter with circular window of 2 meters
f1 <- fastfocal(r, stat = "mean", radius = 2)

# Apply Gaussian-weighted median
f2 <- fastfocal(r, stat = "median", radius = 3, window = "gaussian")
```

---

## Supported Options

| Feature         | Options                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Functions**   | `fastextract()`, `fastfocal()`                                          |
| **Stats**       | `"mean"`, `"min"`, `"max"`, `"sd"`, `"sum"`, `"median"`, `"range"`, `"p25"`, `"p75"` |
| **Windows**     | `"circular"` (default), `"rectangular"`, `"gaussian"`                  |
| **NA Handling** | `na.rm = TRUE` (default) or `FALSE`                                    |
| **Multi-layer** | ✅ Supported                                                           |
| **Progress bar**| ✅ Enabled automatically via `progressr`                               |

---

## Author

Ho Yi Wan  
Associate Professor, University of Florida  
📫 hoyiwan@gmail.com

---

## License

GPL-3
