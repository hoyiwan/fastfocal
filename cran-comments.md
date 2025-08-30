## Test environments
- Local: Windows, R 4.5.1
- win-builder: R-devel
- R-hub: Windows, macOS, Ubuntu runners

## R CMD check results
0 errors | 0 warnings | 0 notes
* New submission.

## Notes
- Vignettes use small examples; heavy benchmark chunks are disabled (eval=FALSE).
- No network access; no files written outside temp directories.
- NA handling matches `terra` (`na.rm`, `na.policy`), FFT auto-switching enabled.
