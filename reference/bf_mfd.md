# Build a brms formula in the congruence–bias parameterisation

Builds the brms non-linear formula for the congruence–bias
parameterisation, in which the native parameters of a
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
or
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
model are functions of the congruence–bias coordinates: the midpoint `M`
of the two expected counts, the congruence `f`, the source bias (through
`delta`, its inverse hyperbolic tangent), and, for `binegbin`, the
overdispersions `kappa`. Each coordinate then takes its own linear
predictor. A predictor such as site in the formula for `logM` alone
multiplies all three rates by a common factor, so the congruence and the
source bias are the same at every site. In the native parameterisation,
the same model needs the predictor in all three rates, with the
coefficients constrained to be equal.

## Usage

``` r
bf_mfd(formula, ...)
```

## Arguments

- formula:

  The response part of the model, written without `~`: `y1 | vint(y2)`
  for
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
  and `y1 | vint(y2, y1_obs)` for
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  and
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).
  The response is passed to the brms formula unchanged.

- ...:

  Two-sided formulae, one for each coordinate, in any order. The
  left-hand side names the coordinate and the right-hand side gives its
  linear predictor, as in `logM ~ 0 + site` or `delta ~ 1`. Every
  coordinate of the chosen set is required, including those with an
  intercept only.

## Value

A `brmsformula` object, as returned by
[`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html).

## Details

The coordinates are related to the native parameters by

\$\$\mu = M f, \quad \lambda_1 = M (1 - f)(1 + \beta), \quad \lambda_2 =
M (1 - f)(1 - \beta), \quad \beta = \tanh \delta,\$\$

where \\\mu\\ is `mu`, the rate of the shared component, and
\\\lambda_1\\ and \\\lambda_2\\ are `lambdaone` and `lambdatwo`, the
rates of the source-exclusive components. Each NB2 (inverse) dispersion
parameter \\\phi\\ is \\1/\kappa^2\\, where \\\kappa\\ is the
corresponding overdispersion. The non-linear parameters are `logM`
(\\\log M\\), `logitf` (\\\operatorname{logit} f\\), `delta`
(\\\delta\\), and the overdispersions. `logM`, `logitf`, and `delta` are
unbounded, so their linear predictors need no link. The overdispersions
must be positive (see Priors).
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
computes the native parameters from the same identities.

**Sets of coordinates.** The set of coordinate formulae supplied
determines the dispersion structure. `logM`, `logitf`, and `delta` are
always required. The three valid sets are:

- `logM`, `logitf`, `delta`, `kappas`, `kappaxone`, and `kappaxtwo`: a
  separate overdispersion for the shared component and for each
  exclusive component, for
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  and
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).

- `logM`, `logitf`, `delta`, `kappas`, and `kappax`: one overdispersion,
  `kappax`, for both exclusive components, so that `shapexone` equals
  `shapextwo`, for
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  and
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).

- `logM`, `logitf`, and `delta`: no overdispersion, for
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md).

Any other set stops with an error. The family passed to
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) must
match the set: a set with overdispersions written for a `bipois` model
declares dpars that the family does not have.

**The generated formula.** For the second set with a site effect on
`logM`,

    bf_mfd(
      y1 | vint(y2),
      logM   ~ 0 + site,
      logitf ~ 1,
      delta  ~ 1,
      kappas ~ 1,
      kappax ~ 1
    )

returns the same formula as

    bf(y1 | vint(y2) ~ 1, nl = TRUE) +
      nlf(mu        ~ logM + log_inv_logit(logitf)) +
      nlf(lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)) +
      nlf(lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)) +
      nlf(shapes    ~ -2 * log(kappas)) +
      nlf(shapexone ~ -2 * log(kappax)) +
      nlf(shapextwo ~ -2 * log(kappax)) +
      lf(logM ~ 0 + site, logitf ~ 1, delta ~ 1, kappas ~ 1, kappax ~ 1)

The first set replaces the two `kappax` lines with
`shapexone ~ -2 * log(kappaxone)` and `shapextwo ~ -2 * log(kappaxtwo)`;
the third set omits all three dispersion lines.

**Priors.** Every coordinate is a non-linear parameter, to which brms
gives a flat, improper prior unless one is set. We recommend a prior on
each, written with `class = "b"` and `nlpar` set to the name of the
coordinate. The expression \\-2 \log \kappa\\ is undefined for \\\kappa
\le 0\\, so a prior on an overdispersion with an intercept-only formula
needs `lb = 0`, as in
`prior(normal(0, 1), class = "b", nlpar = "kappas", lb = 0)`. The worked
example in
[`vignette("bicountbrms")`](https://anhsmith.github.io/bicountbrms/articles/bicountbrms.md)
gives priors for every coordinate.

`bf_mfd()` adds no family, Stan function, or post-processing method. The
fitted model is a `bipois` or `binegbin` fit whose native parameters are
non-linear functions of the coordinates, and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
and [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
dispatch to the family methods as for any other fit.

## See also

[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md),
[`dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/dpars_to_mfd.md),
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)

## Examples

``` r
# A site effect on the midpoint alone, with one overdispersion for both
# exclusive components
bf_mfd(
  y1 | vint(y2),
  logM   ~ 0 + site,
  logitf ~ 1,
  delta  ~ 1,
  kappas ~ 1,
  kappax ~ 1
)
#> y1 | vint(y2) ~ 1 
#> mu ~ logM + log_inv_logit(logitf)
#> lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)
#> lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)
#> shapes ~ -2 * log(kappas)
#> shapexone ~ -2 * log(kappax)
#> shapextwo ~ -2 * log(kappax)
#> logM ~ 0 + site
#> logitf ~ 1
#> delta ~ 1
#> kappas ~ 1
#> kappax ~ 1

# A partially paired design with Poisson components
bf_mfd(
  y1 | vint(y2, y1_obs),
  logM   ~ 0 + site,
  logitf ~ 1,
  delta  ~ 1
)
#> y1 | vint(y2, y1_obs) ~ 1 
#> mu ~ logM + log_inv_logit(logitf)
#> lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)
#> lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)
#> logM ~ 0 + site
#> logitf ~ 1
#> delta ~ 1
```
