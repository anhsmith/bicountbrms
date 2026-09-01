# Migration and errata

## Package history

These families were first released as `skellambrms` (0.1.0–0.5.0) and
then as `pairedcountbrms` (0.6.0–0.8.0). That second package held two
unrelated suites: families modelling a count pair jointly, and families
modelling its difference $`d = y_1 - y_2`$ directly. The two shared no
code. Release 0.9.0 separated them. The joint families became this
package; the difference families returned to
[`skellambrms`](https://github.com/anhsmith/skellambrms).

The separation renamed nothing, so
[`library(pairedcountbrms)`](https://rdrr.io/r/base/library.html)
becomes
[`library(bicountbrms)`](https://github.com/anhsmith/bicountbrms) or
[`library(skellambrms)`](https://rdrr.io/r/base/library.html) according
to which suite was in use. `github.com/anhsmith/pairedcountbrms`
redirects here.

## Changes by release

**0.8.0** replaced the single excess dispersion `shapex` with the
per-margin pair `shapexone`/`shapextwo`, for the partially observed
family alone.

**0.9.0** renamed `binegbin_joint()` to `binegbin_cens()`, on the
grounds that every family in the package modelled the pair jointly and
`_joint` therefore distinguished nothing. The same release added
`bipois_cens()`.

**0.10.0** replaces four families with two. Each of `bipois` and
`binegbin` now has two constructors: a plain one for a fully paired
design, and a `_partialobs()` one for a design in which the first count
is missing on some rows. Both constructors of a pair return the same
`custom_family` name, so one Stan function and one set of
post-processing methods serve both. The same release gives
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
the six distributional parameters previously held by `binegbin_cens()`
alone.

## Substitutions for existing code

| up to 0.9.1 | from 0.10.0 |
|----|----|
| `bipois_cens()` | [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md) |
| `binegbin_cens()`, `binegbin_joint()` | [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md) |
| `bipois_cens_stanvars()` | [`bipois_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md) |
| `binegbin_cens_stanvars()`, `binegbin_joint_stanvars()` | [`binegbin_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md) |
| [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md) with `shapex ~ 1` | [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md) with `shapexone ~ 1, shapextwo ~ 1` |

The old names are removed outright. There is no deprecation layer, and
nothing is left to remove at a later release.

A formula written against the five-dpar
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
needs its dispersion term replaced. Either estimate the two separately:

``` r

bf(y1 | vint(y2) ~ 1, mu ~ 1,
   nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx),
   lamx ~ 1, shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE)
```

or reproduce the pre-0.10.0 model exactly, by routing both through one
non-linear parameter:

``` r

bf(y1 | vint(y2) ~ 1, mu ~ 1,
   nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx),
   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
   lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE)
```

The second form is term for term the likelihood the single `shapex`
gave, and a package test verifies that the generated Stan assigns both
distributional parameters from one non-linear parameter.

Under the tied form, the dispersion prior is written

``` r

prior(normal(0, 1.5), class = "b", nlpar = "shapexx")
```

rather than the pre-0.10.0 `class = "Intercept", dpar = "shapex"`. Both
fields change. A prior written the old way names no parameter in the new
model and is dropped without a warning, which leaves the dispersion
improper. See [Choosing
priors](https://anhsmith.github.io/bicountbrms/articles/choosing-priors.md).

## `_cens` named the wrong mechanism

brms already uses `cens()` as an addition term meaning a value known to
lie in a set — `left`, `right`, `interval` — and supplies no code for a
value that was not observed at all. On the rows these constructors
admit, the first count was not recorded, and the likelihood marginalises
over its whole support rather than over a bounded set. A brms user
meeting `_cens` had every reason to map one mechanism onto the other.

`partialobs` names what is partial, which is the *pair*, and which is a
property of the data rather than of the model. `unobs` and `unmatched`
name only the minority class of rows; `partial` alone collides with
partial likelihood and partial pooling.

## Erratum: `posterior_epred_binegbin()` before 0.9.0

Releases before 0.9.0 computed
[`posterior_epred_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
using the *marginal* shared fraction $`\mu/(\mu + \lambda_2)`$ in place
of the *conditional* fraction. That substitution is `bipois`’s answer,
exact only in the Poisson limit, and it made the expectation disagree
with the family’s own
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
draws by more than Monte Carlo error.

Over a grid of plausible rates and dispersions, with $`y_2`$ within one
standard deviation of its mean, the substituted value was in error by
more than 5% in about half the settings and by more than 20% in a third
of them. The direction of the error depends on which component is the
more dispersed.

[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
were not affected at any release. Any number computed with an earlier
[`posterior_epred_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
should be recomputed.

## Erratum: the dispersion routing in `posterior_predict()`, before 0.9.1

Release 0.9.1 added `tests/testthat/test-cens-predict.R`, since
superseded by `test-partialobs-predict.R`, after a check of the test
suite’s discrimination. In `posterior_predict_binegbin_cens()`, the
conditional split of $`N_{\text{shared}} \mid y_2`$ is weighted by
`shapextwo`, while the fresh source-specific count added on top takes
`shapexone`. Exchanging those two passed 412 assertions across six test
files without a single failure, because every prediction test then in
the suite built its fixture from a single five-dpar `shapex`, which both
name vectors resolve to identically.

No released version computed the exchanged quantity. The defect was in
the test suite’s power to detect an error, not in the shipped code.

## References
