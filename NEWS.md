# fastfocal 0.1.0

* Initial CRAN release.

## Features
- Fast focal / moving-window operations for `SpatRaster` with **auto engine**:
  C++ via `terra` for small windows, **FFT** for large.
- **NA handling matches `terra`**:
  - `na.rm = FALSE` uses full-box gating (any NA in the kr×kc box → NA).
  - `na.rm = TRUE` averages over available cells in the kernel support.
  - `na.policy = "omit"` or `"all"` supported; center NA behavior aligned with `terra`.
- **Kernel support**: `"circle"`, `"rectangle"`, `"gaussian"`, `"pareto"`, `"idw"`,
  `"exponential"`, `"triangular"`, `"cosine"`, `"logistic"`, `"cauchy"`,
  `"quartic"`, `"epanechnikov"`, or pass a **numeric matrix** as the kernel.
- **FFT padding (`pad = "auto"`)** to avoid prime-factor slowdowns; uses next-fast length.
- Smarter `engine="auto"` heuristic; accepts multi-layer rasters.

## Extraction
- `fastextract()` for points/polygons with optional buffer radii; returns tidy summaries.

## Docs
- Two vignettes (overview + benchmarks). Heavy benchmark chunks are disabled for CRAN.
