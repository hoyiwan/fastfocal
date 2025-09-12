## fastfocal 0.1.2 — CRAN submission

### Test environments
- Local: Windows 11 x64, R 4.5.1 (ucrt)
- win-builder: R-devel (Windows Server 2022)
- R-hub GitHub Actions: Windows, macOS (Intel & arm64), Ubuntu release

### R CMD check results
- Local: 0 errors | 0 warnings | 0 notes
- win-builder (R-devel): 0 errors | 0 warnings | 1 note  
  *NOTE:* “New submission” and “Possibly misspelled words in DESCRIPTION: ‘Hijmans’, ‘terra’.”  
  These are a surname (Hijmans) and an R package name (‘terra’), intentionally spelled.
- R-hub: 0 errors | 0 warnings | 0 notes (all selected platforms)

### Changes and policy compliance
- Software/package names are single-quoted in the Description as per CRAN policy (e.g., ‘terra’). The Title avoids quotes.
- References are cited in vignettes and documentation, not in DESCRIPTION.
- Examples and vignettes:
  - Small and fast; heavy benchmark code disabled for CRAN checks.
  - Save and restore `par()` to avoid side effects.
  - No network access; no files written outside temp directories.
- Functionality notes:
  - NA handling matches ‘terra’ (`na.rm`, `na.policy = "omit" | "all"`).
  - `engine = "auto"` selects ‘C++’ or ‘FFT’ depending on scale; optional FFT padding to avoid prime-factor slowdowns.
  - Matrix kernels accepted, matching `terra::focalMat()` where needed.
- Packaging:
  - UTF-8 encoding; no non-ASCII in code.
  - `.Rbuildignore` excludes `.github/`, `cran-comments.md`, and `CITATION.cff`.

### Reverse dependencies
- None.

### Additional notes
- CRAN URL (https://CRAN.R-project.org/package=fastfocal) currently returns 404 as expected; it will resolve once the package is published.
