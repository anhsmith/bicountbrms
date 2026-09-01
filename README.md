# bicountbrms <img src="man/figures/logo.png" align="right" height="139" alt="bicountbrms logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Custom [brms](https://paulbuerkner.com/brms/) families for Bayesian modelling of
**pairs of counts** that are positively correlated. Applications include counts
recorded by two sources that measure the same quantity, such as two observers, two
instruments, or two reporting channels, where the question is how much, 
and how systematically, the two sources disagree.

The pair $(y_1, y_2)$ is modelled jointly using **trivariate reduction** 
with Poisson (Holgate 1964; Karlis and Ntzoufras 2003) 
or negative-binomial (Kirkpatrick and Neale 2016; Kirkpatrick 2022) components, 
using `bipois()` or `binegbin()`, respectively.
Correlation between the two counts is induced via an unobserved shared latent component, 
with an additional component each that only one source recorded. 

The native parameters
— means and (for negative binomial) dispersions for each of the three
latent components — can be reparameterised into four interpretable
coordinates: the overall count rate, the congruence between the two
sources, the bias toward one of them, and (for negative binomial)
the excess dispersion of each source.
The package also allows the model to be fit when one source is only 
partially observed (the `*_partialobs` variants); posterior predictions 
can then be made for the missing values of that source 
based on the observed values of the other. 


## The core functions

|  | both counts recorded | first count sometimes unrecorded |
|---|---|---|
| **equidispersed** (Poisson latent counts) | `bipois()` | `bipois_partialobs()` |
| **overdispersed** (negative-binomial latent counts) | `binegbin()` | `binegbin_partialobs()` |

Each row of the table above is one family (Poisson or negative binomial). 
The two constructors for each family return the same 
likelihood and the same post-processing methods, and differ only in whether the 
$y_1$ observation indicator (`y1_obs`) enters the likelihood as data or as a constant. 
The indicator is not supplied in cases where both counts are fully observed.

For a row where both counts are observed, the likelihood term is the
full joint distribution of $(y_1, y_2)$: a function of the shared rate
and (for negative binomial) dispersion (`mu`, `shapes`), each source's own
rate and (for negative binomial) dispersion (`lambdaone`, `shapexone`,
`lambdatwo`, `shapextwo`), and, through
`lambdaone` and `lambdatwo`, the bias between the two sources.

For a row where only $y_2$ is observed, the likelihood term is instead
the marginal probability of $y_2$ alone, obtained by summing that joint
distribution over every value the first count could have taken. 
It therefore only informs the parameters common to both counts (`mu`,
`shapes`, `lambdatwo`, `shapextwo`) and not the parameters that pertain 
only to $y_1$ (`lambdaone`, `shapexone`).


## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/bicountbrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Either `rstan` or `cmdstanr` works as the `brms` backend.

## Minimal fits

The second count is provided as supplementary integer data via the `vint()` function, and
the family's Stan function is specified via the `stanvars` argument. Both are required.

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

Where the first count was not recorded on some rows, each row supplies a second
integer through `vint()`, in the order `vint(y2, y1_obs)`, and the family and its
stanvars are the `_partialobs` constructors:

```r
fit_po <- brm(
  bf(y1 | vint(y2, y1_obs) ~ 1,
     lambdaone ~ 1, lambdatwo ~ 1,
     shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
  family   = binegbin_partialobs(),
  stanvars = binegbin_partialobs_stanvars(),
  data     = dat
)
```

`y1_obs` is `1` where the first count was recorded and `0` where it was not. On a
row with `y1_obs == 0`, `y1` may hold any non-negative integer — `0` is the
conventional placeholder — because that row's likelihood term does not depend on
it. Do not use `NA`: brms drops such rows before fitting, taking their observed
$y_2$ with them. `posterior_predict()` and `posterior_epred()` then return a draw
and $E[y_1 \mid y_2]$ for every row, including the rows where $y_1$ was never
recorded.

`brms` leaves every distributional parameter of a custom family except `mu`
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
| `y1_obs` (via `vint()`) | — | the `_partialobs()` indicator: `1` where $y_1$ was recorded and `0` otherwise |

`brms::custom_family()` requires one distributional parameter to be named `mu`,
and rejects any name ending in a digit. Here `mu` is the rate of the shared
component rather than the mean of either count, and the two source-specific
rates are spelled out (`one` or `two`).

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

**Modelling the differences directly.** Where only the disagreement
$d = y_1 - y_2 \in \mathbb{Z}$ is of interest, with optional truncation via
`resp_trunc()`, the companion package
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
Ntzoufras (2003). Kirkpatrick and Neale (2016) and Kirkpatrick (2022) use the same construction
with negative-binomial components. Full references are given in [The families and
their parameters][families].

[started]: https://anhsmith.github.io/bicountbrms/articles/bicountbrms.html
[families]: https://anhsmith.github.io/bicountbrms/articles/families-and-parameters.html
[priors]: https://anhsmith.github.io/bicountbrms/articles/choosing-priors.html
[partial]: https://anhsmith.github.io/bicountbrms/articles/partially-observed-fit.html
[anatomy]: https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.html
[migration]: https://anhsmith.github.io/bicountbrms/articles/migration-and-errata.html
