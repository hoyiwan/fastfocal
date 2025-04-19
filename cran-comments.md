## Test environments
* local Windows 11, R 4.4.3
* devtools::check(), devtools::check_built() all passed
* Manual vignette rendering used (via rmarkdown::render) due to known issues with Quarto detection

## R CMD check results
There were no ERRORs or WARNINGs.
There was 1 NOTE:
* "unable to verify current time" — this is a known and harmless time zone check on Windows.

## Vignettes
* Vignettes were rendered manually into inst/doc/
* `devtools::build(vignettes = FALSE)` was used to avoid Quarto-related check failures.
