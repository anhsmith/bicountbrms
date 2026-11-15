# bicountbrms <img src="man/figures/logo.png" align="right" height="139" alt="bicountbrms logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22239120.svg)](https://doi.org/10.5281/zenodo.22239120)
<!-- badges: end -->

Custom [brms](https://paulbuerkner.com/brms/) families for Bayesian modelling of
**pairs of counts** that are positively correlated. Applications include counts
recorded by two sources that measure the same quantity, such as two observers, two
instruments, or two reporting channels.

The pair $(y_1, y_2)$ is modelled by **trivariate reduction**, with Poisson
components in `bipois()` (Holgate 1964; Karlis and Ntzoufras 2003) or
negative-binomial components in `binegbin()` (Kirkpatrick and Neale 2016;
Kirkpatrick 2022). Each count is the sum of a latent count recorded by both
sources, which induces the correlation, and a latent count recorded by that
source alone.

The native parameters of the model are the rates and (for negative binomial)
dispersions of the three latent components. The model can also be fitted using
alternative [congruence–bias coordinates][reparameterisation], each of which
directly represents a quantity of interest when comparing two sources: the
midpoint of the two expected counts, the congruence between the two sources,
the source bias, and (for negative binomial) the overdispersion of each
component. `bf_mfd()` writes the brms formula in these coordinates, so priors
and predictors act directly on the quantities of interest, and `summary(fit)`
reports the posterior of each coordinate. The native
parameterisation can be fitted with an ordinary brms formula, in which a
predictor can instead act on a single latent component; see *Models with
predictors* in [Getting started][started].

The `_partialobs()` constructors fit the model when one source has missing
values, and the fitted model can impute the missing counts of that source
conditional on the recorded counts from the other source.


## The core functions

|  | both counts recorded | source 1 has missing values |
|---|---|---|
| **equidispersed** (Poisson latent counts) | `bipois()` | `bipois_partialobs()` |
| **overdispersed** (negative-binomial latent counts) | `binegbin()` | `binegbin_partialobs()` |

Each row of the table is one family. The two constructors of a family share the
likelihood and the post-processing methods, and the `_partialobs()` constructor
also takes the observation indicator `y1_obs`.


## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/bicountbrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Either `rstan` or `cmdstanr` works as the `brms` backend.

## Minimal fits

The models below are written in the [congruence–bias coordinates][reparameterisation]
with `bf_mfd()`, which takes one formula for each coordinate (`logM`, `logitf`,
`delta`, `kappas`, `kappaxone`, and `kappaxtwo`; see [Getting started][started]
for their definitions).
$y_2$ is passed to the family as supplementary integer data through
`vint()`, and the Stan function for the family through the `stanvars`
argument. Both are required. brms gives every coordinate a flat, improper prior
unless one is set, and [Getting started][started] gives a prior for each
coordinate.

### Paired designs

```r
library(brms)
library(bicountbrms)

fit <- brm(
  bf_mfd(
    y1 | vint(y2),
    logM      ~ 1,
    logitf    ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  ),
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  data     = dat,
  prior    = priors  # one prior for each coordinate; see Getting started
)
```

### Partially paired designs

Where source 1 has missing values, the data need an observation indicator
`y1_obs`, `1` where $y_1$ was recorded and `0` where it was not. **`y1` must not
be `NA`**, because brms drops any row with a missing response with only a
warning. On an unpaired row, the likelihood term does not depend on `y1`, so
`y1` can be any integer. We suggest `-1`, because it cannot be confused with a
recorded count:

```r
dat$y1_obs <- as.integer(!is.na(dat$y1))
dat$y1[is.na(dat$y1)] <- -1L
```

The indicator is passed as the second argument of `vint()`, and the family and
its stanvars are the `_partialobs` constructors:

```r
fit_po <- brm(
  bf_mfd(
    y1 | vint(y2, y1_obs),
    logM      ~ 1,
    logitf    ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  ),
  family   = binegbin_partialobs(),
  stanvars = binegbin_partialobs_stanvars(),
  data     = dat,
  prior    = priors
)
```

To impute $y_1$ on the rows where it was not recorded, use
`posterior_predict_y1_given_y2()`; see [Prediction and imputation][prediction].

## Documentation

- [Getting started][started] describes one `binegbin()` analysis from start to
  finish: simulating data from known parameters, fitting the model in the
  congruence–bias parameterisation, checking recovery of the generating values,
  predicting both counts, and comparing the fit with a Poisson model.
- [Partially paired designs][partial] fits a model in the congruence–bias
  parameterisation to simulated data in which source 1 has missing values, and
  compares the imputed values of $y_1$, given $y_2$ and marginal, with values
  withheld from the model.
- [Prediction and imputation][prediction] describes the prediction
  functions and the imputation of $y_1$ on unpaired rows.
- [The congruence–bias parameterisation][reparameterisation] has an interactive
  widget that simulates pairs of counts from one model and shows the
  congruence–bias and native parameterisations side by side, and gives the
  equations that convert between the coordinates and the native parameters.
- The [changelog][news] records what each release changed.

A companion package,
[`skellambrms`](https://github.com/anhsmith/skellambrms), supplies Skellam,
discrete Laplace, and discrete normal families that model $d = y_1 - y_2$ directly, 
with optional truncation via `resp_trunc()`.

## Funding

This package was developed by Sea Through Science for a fisheries
electronic-monitoring project led by Johanna Pierre of JPEC Consulting, with
support from The Nature Conservancy, the Belize Blue Bonds for Ocean
Conservation, and the Belize Fund for a Sustainable Future.

## References

Holgate, P. (1964). Estimation for the bivariate Poisson distribution.
*Biometrika* **51**(1–2), 241–245. <https://doi.org/10.2307/2334210>

Karlis, D. and Ntzoufras, I. (2003). Analysis of sports data by using bivariate
Poisson models. *Journal of the Royal Statistical Society: Series D (The
Statistician)* **52**(3), 381–393. <https://doi.org/10.1111/1467-9884.00366>

Kirkpatrick, R. M. (2022). *RMKdiscrete: Sundry discrete probability
distributions*. R package version 0.2.
<https://CRAN.R-project.org/package=RMKdiscrete>

Kirkpatrick, R. M. and Neale, M. C. (2016). Applying multivariate discrete
distributions to genetically informative count data. *Behavior Genetics*
**46**(2), 252–268. <https://doi.org/10.1007/s10519-015-9757-z>

[started]: https://anhsmith.github.io/bicountbrms/articles/bicountbrms.html
[partial]: https://anhsmith.github.io/bicountbrms/articles/partially-paired-designs.html
[prediction]: https://anhsmith.github.io/bicountbrms/articles/prediction.html
[reparameterisation]: https://anhsmith.github.io/bicountbrms/articles/reparameterisation.html
[news]: https://anhsmith.github.io/bicountbrms/news/index.html
