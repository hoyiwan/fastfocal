## Test environments
- Local: R --as-cran (0 errors | 0 warnings | 0 notes)
- GitHub Actions: macOS, Windows, Ubuntu (R-devel, R-release, oldrel)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Resubmission
- Address CRAN request to reset graphics state in vignettes:
  - Replaced on.exit() in top-level chunks with explicit save/restore:
    * `oldpar <- par(no.readonly = TRUE)` … `par(oldpar)`
    * `layout(1)` after using `layout()`
- No API or functional changes to exported R code.
