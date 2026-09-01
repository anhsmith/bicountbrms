# bicountbrms <img src="man/figures/logo.png" align="right" height="139" alt="bicountbrms logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Custom [brms](https://paulbuerkner.com/brms/) families for a **pair of counts**
recorded by two sources that measure the same quantity — two observers, two
instruments, two reporting channels — where the question is how much, and how
systematically, the two disagree.

The pair $(y_1, y_2)$ is modelled jointly rather than reduced to its
difference, using **trivariate reduction** (Holgate 1964): the two counts share
an unobserved latent component, which induces their correlation, plus a
component each that only one source recorded. One model then estimates the
overall level of counting, the congruence between the two sources, the bias
toward one of them, and the excess dispersion of each. A fitted model also
predicts either count from the other, including on rows where one of them was
never recorded.

## The four constructors

|  | both counts recorded | first count sometimes unrecorded |
|---|---|---|
| **equidispersed** (Poisson latent counts) | `bipois()` | `bipois_partialobs()` |
| **overdispersed** (negative-binomial latent counts) | `binegbin()` | `binegbin_partialobs()` |

Each row of that table is one family. Both of its constructors return the same
likelihood and the same post-processing methods, and differ only in whether the
data include an observation flag. A fully paired design never supplies one.

Rows on which the first count was not recorded are scored by the second count's
marginal, taken from the same joint model, so they still inform the shared
component and the second source rather than being dropped. What they cannot
inform is the first source's own rate and dispersion, or the bias between the
two sources: those appear only where both counts were recorded, so the matched
subset is the effective sample size for them. That is worth knowing before the
data are collected. This mechanism is unrelated to brms's `cens()` addition
term, which describes a value known to lie in a set.

## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/bicountbrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Either rstan or cmdstanr works as the brms backend.

## A minimal fit

The second count travels as supplementary integer data through `vint()`, and
the family's Stan function arrives through `stanvars`. Both are required.

```r
library(brms)
library(bicountbrms)

fit <- brm(
  bf(y1 | vint(y2) ~ 1,
     lambdaone ~ 1, lambdatwo ~ 1,
     shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  data     = dat
)
```

brms leaves every distributional parameter of a custom family except `mu`
without a prior. Set them; see [Choosing priors][priors].

## Notation

| Code | Math | Meaning |
|---|---|---|
| `y1` (the response) | $y_1$ | count from source 1 |
| `y2` (via `vint()`) | $y_2$ | count from source 2 |
| — | $N_{\text{shared}}$ | latent count both sources recorded |
| — | $N_1$, $N_2$ | latent counts only one source recorded |
| `mu` | $\mu$ | rate of the shared component (**not** a response mean) |
| `lambdaone`, `lambdatwo` | $\lambda_1$, $\lambda_2$ | rates of the two source-specific components |
| `shapes` | $\phi_{\text{s}}$ | NB2 dispersion of the shared component |
| `shapexone`, `shapextwo` | $\phi_{x1}$, $\phi_{x2}$ | NB2 dispersions of the two source-specific components |
| `y1_obs` (via `vint()`) | — | the `_partialobs()` flag: was $y_1$ recorded? |

`brms::custom_family()` requires one distributional parameter to be named `mu`,
and rejects any name ending in a digit. Here `mu` is the rate of the shared
component rather than the mean of either count, and the two source-specific
rates are spelled out. Source 1 is whichever count is placed on the left-hand
side of the formula; the labelling means nothing further.

## Documentation

- [Getting started][started] — one fit, end to end: simulate, fit, check,
  predict.
- [The families and their parameters][families] — the construction, the
  moments, and a fitted demonstration that the two excess dispersions can be
  told apart when they differ.
- [Choosing priors][priors] — what brms leaves improper, and which of the
  `class`, `dpar` and `nlpar` slots each prior belongs in.
- [A worked partially observed fit][partial] — imputing a count that was never
  recorded, scored against values withheld from the model.
- [The anatomy of a paired count][anatomy] — overall level, congruence and
  source bias as coordinates, with shrinkage priors.
- [Migration and errata][migration] — what each release changed, and two
  corrected results.

**Modelling the difference instead.** Where only the disagreement
$d = y_1 - y_2 \in \mathbb{Z}$ is of interest, or where truncation via
`resp_trunc()` is required, the companion package
[`skellambrms`](https://github.com/anhsmith/skellambrms) supplies Skellam,
discrete Laplace and discrete normal families that model $d$ directly. Both
packages were formerly distributed together as `pairedcountbrms`; see
[Migration and errata][migration].

## Citation

```r
citation("bicountbrms")
```

Holgate (1964) states the trivariate-reduction representation and gives
maximum-likelihood estimation for it; the regression form is due to Karlis and
Ntzoufras (2003). Kirkpatrick and Neale (2016, 2022) use the same construction
with negative-binomial components. Full references are given in [The families and
their parameters][families].

[started]: https://anhsmith.github.io/bicountbrms/articles/bicountbrms.html
[families]: https://anhsmith.github.io/bicountbrms/articles/families-and-parameters.html
[priors]: https://anhsmith.github.io/bicountbrms/articles/choosing-priors.html
[partial]: https://anhsmith.github.io/bicountbrms/articles/partially-observed-fit.html
[anatomy]: https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.html
[migration]: https://anhsmith.github.io/bicountbrms/articles/migration-and-errata.html
