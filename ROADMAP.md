# Future work

Proposed additions to bicountbrms, not yet scheduled. Each entry states
the change and the design questions it raises. Move an entry to
`NEWS.md` when it is released.

## Native and (M, f, δ) parameterisations throughout the documentation

The get-started vignette would introduce both parameterisations: the
native rates and dispersions (`mu`, `lambdaone`, `lambdatwo`, `shapes`,
`shapexone`, `shapextwo`), and the midpoint $`M`$, congruence $`f`$ and
source bias $`\delta`$ converted by
[`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md)
and
[`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md).
The remaining articles would then give each formula, prior and summary
in both parameterisations, as tabs (`## Heading {.tabset}` in R
Markdown, which pkgdown renders).

- Decide whether each Mfd tab fits its own model, which adds one fit per
  tab to the precompile, or converts the posterior draws of the native
  fit.
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

- The observation flags would become two `vint()` integers, so the plain
  constructors would pass two literals (`"1"`, `"1"`) rather than one.
  Every constructor changes its `vars`, and
  `tests/testthat/test-stancode-shape.R` must be updated for all four.
- `.y1_obs_at()` in `R/utils.R` generalises to a lookup for either flag,
  and `posterior_epred_*` needs $`E[y_2 \mid y_1]`$ for rows where only
  the first count was recorded.
- The split between the shared and excess components is estimated mainly
  from the paired rows, as it is now, so the number of paired rows
  remains the effective sample size for the bias and the congruence.
