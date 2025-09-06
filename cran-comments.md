## fastfocal 0.1.1 — CRAN submission

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
- Software/package names are single-quoted in Title/Description (‘FFT’, ‘C++’, ‘terra’).
- Method references included in DESCRIPTION:
  - Cooley & Tukey (1965) <doi:10.1090/S0025-5718-1965-0178586-1>
  - Hijmans, R. J. (2024). ‘terra’ R package. <https://CRAN.R-project.org/package=terra>
- Examples and vignettes:
  - Use small, fast, executable snippets.
  - Save and restore `par()` to avoid side effects.
  - No network access; no files written outside temp directories.
- Functionality notes:
  - NA handling matches ‘terra’ (`na.rm`, `na.policy = "omit" | "all"`).
  - `engine = "auto"` selects ‘C++’ or ‘FFT’ based on scale; optional padding to FFT-friendly sizes.
  - Accepts matrix kernels to match `terra::focalMat()` exactly when needed.
- Packaging:
  - UTF-8 encoding; no non-ASCII in code.
  - `.Rbuildignore` includes `.github/` and `cran-comments.md`.

### Reverse dependencies
- None.
