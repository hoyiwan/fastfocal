# fastfocal 0.1.3

This release represents the version used to generate all benchmark results and figures reported in the forthcoming manuscript.

Improvements
- Performance optimizations
- Refined core routines for faster runtime in typical benchmark scenarios.

# fastfocal 0.1.2

* Fix: NA-semantic parity with `terra::focal()` under `na.rm=TRUE/FALSE` and
  `na.policy="omit"`; removed speckle artifacts at NA boundaries.
* Fix: Correct matrix orientation and final transpose in FFT path.
* Perf: Auto-padding to next 5-smooth lengths (`2*3*5` factors) to avoid large-prime
  slowdowns in FFTs.
* Rebuilt vignettes and refreshed documentation.
* Fixed Rd markup (e.g., replaced `{terra}` with `\pkg{terra}`).
* Removed non-ASCII characters from R sources and roxygen; math now uses Rd math (e.g., `\eqn{\sigma}`).
* Vignettes now restore user graphics state:
  - Added `oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)` in chunks that modify `par(...)`.
  - Added `on.exit(layout(1), add = TRUE)` in chunks that use `layout(...)`.
* No changes to `options()` or working directory in examples/vignettes.

# fastfocal 0.1.1

* First public release (submitted to CRAN).
* Added Zenodo DOI and updated citation information.
* Updated vignettes (`index.Rmd` and `benchmark.Rmd`) with correct version numbers and DOI.
* Clarified that the FFT backend supports only `sum` and `mean`; other statistics handled by the `terra` backend.
* README and DESCRIPTION updated for consistency (license MIT, DOI, installation instructions).
* Added `inst/CITATION` file so that `citation("fastfocal")` prints the DOI and pkgdown URL.

# fastfocal 0.1.0

* Internal development version (not submitted to CRAN).

# Citation
- If using fastfocal in academic work, cite this version as:

Wan, H.Y. (2025). fastfocal v0.1.3. GitHub. https://github.com/yourrepo/fastfocal/releases/tag/v0.1.3
