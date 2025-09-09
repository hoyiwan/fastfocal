# fastfocal 0.1.2

* Fix: NA-semantic parity with `terra::focal()` under `na.rm=TRUE/FALSE` and
  `na.policy="omit"`; removed speckle artifacts at NA boundaries.
* Fix: Correct matrix orientation and final transpose in FFT path.
* Perf: Auto-padding to next 5-smooth lengths (`2*3*5` factors) to avoid large-prime
  slowdowns in FFTs.

# fastfocal 0.1.1

* First public release (submitted to CRAN).
* Added Zenodo DOI and updated citation information.
* Updated vignettes (`index.Rmd` and `benchmark.Rmd`) with correct version numbers and DOI.
* Clarified that the FFT backend supports only `sum` and `mean`; other statistics handled by the `terra` backend.
* README and DESCRIPTION updated for consistency (license MIT, DOI, installation instructions).
* Added `inst/CITATION` file so that `citation("fastfocal")` prints the DOI and pkgdown URL.

# fastfocal 0.1.0

* Internal development version (not submitted to CRAN).
