# Future work

Proposed additions to bicountbrms, not yet scheduled. Each entry states
the change and the design questions it raises. Move an entry to
`NEWS.md` when it is released.

## Both parameterisations in the articles

Getting started introduces both parameterisations. The articles would
give each formula, prior and summary in both the congruence–bias
coordinates and the native parameters, as tabs (`## Heading {.tabset}`
in R Markdown, which pkgdown renders).

- Decide whether each congruence–bias tab fits its own model, which adds
  one fit per tab to the precompile, or converts the posterior draws of
  the native fit.
- Priors written on one scale induce a different density on the other,
  so the two tabs are not the same model unless the prior is converted
  as well.

## Zero inflation

Excess zeros in one or both observed counts, beyond those implied by the
negative-binomial components.

- Decide where the zero process acts: on each observed count separately,
  on the shared component, or on the pair (a structural zero for both
  counts). These are three different models, and only the last one
  leaves the shared component interpretable as a count recorded by both
  sources.
- A zero-inflation probability is a new distributional parameter, which
  changes the lpmf signature for every constructor that gains it.

## Both sources partially observed

At present only the first count may be unrecorded. Allowing the second
to be unrecorded as well adds a third kind of row, in which only the
first count is observed, scored by the marginal of the first count.

- The observation indicators would become two `vint()` integers, so the
  plain constructors would pass two literals (`"1"`, `"1"`) rather than
  one. Every constructor changes its `vars`, and
  `tests/testthat/test-stancode-shape.R` must be updated for all four.
- `.y1_obs_at()` in `R/utils.R` generalises to a lookup for either
  indicator, and new functions `posterior_predict_y2_given_y1()` and
  `posterior_epred_y2_given_y1()` would impute y_2 on rows where only
  y_1 was recorded.
- The split between the shared and exclusive components is estimated
  mainly from the paired rows, as it is now, so the number of paired
  rows remains the effective sample size for the bias and the
  congruence.

## Citation section in the README

The Citation section was removed from `README.md` on 14 November 2026,
before the paper describing the package was submitted to *Methods in
Ecology and Evolution*, because the reviewers receive an anonymised copy
of the repository. Restore the Citation section once the paper is
published.

- Point the restored section and `inst/CITATION` to the paper, and keep
  the `citation("bicountbrms")` call and the concept DOI
  (10.5281/zenodo.22239120).
- Add Kirkpatrick (2022) to the `citHeader` in `inst/CITATION`, which at
  present cites Holgate (1964), Karlis and Ntzoufras (2003), and
  Kirkpatrick and Neale (2016).
