## Submission

This is a new submission.

## Test environments

* Local: Windows 11, R 4.6.0 (2026-04-24 ucrt), `R CMD check --as-cran`.
* win-builder: R-devel (2026-09-30 r90605 ucrt): 1 NOTE, described below.
* mac builder: macOS 26.6 (arm64), R 4.6.1 Patched (2026-07-27 r90311): Status OK.
* GitHub Actions: ubuntu-latest, R release.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* Possibly misspelled words in DESCRIPTION: "trivariate" and
  "reparameterisation" are correctly spelled technical terms (trivariate
  reduction is the construction the package implements). The function names
  flagged on win-builder are now written with parentheses.
* Suggests or Enhances not in mainstream repositories: cmdstanr.

  cmdstanr is an optional alternative to rstan as the 'brms' backend. It is
  distributed from <https://stan-dev.r-universe.dev>, which is declared in
  `Additional_repositories`. Neither the package code nor the tests call it;
  it is used only by the vignette sources, which are excluded from the build.
  The package works with rstan alone, which 'brms' imports.

## Check time and Stan compilation

Fitting a model with these families compiles a Stan program, which takes
minutes per model. The tests that compile Stan run only when `NOT_CRAN=true`:
those that call `brm()` are guarded by `skip_on_cran()`, and the file-level
`rstan::stan_model()` blocks by an explicit test of `NOT_CRAN`. On CRAN, no
Stan program is compiled. The local build and `--as-cran` check took 95 s, of
which the tests took 37 s.

The vignette is precompiled: `vignettes/bicountbrms.Rmd` is knitted in
advance from a source file excluded by `.Rbuildignore`, so building it fits no
models.

## Reverse dependencies

There are no reverse dependencies.
