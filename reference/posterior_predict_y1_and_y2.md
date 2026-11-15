# Draws and expectations of both counts from the joint model

`posterior_predict_y1_and_y2()` draws the pair `(y1, y2)` from the joint
model of the fit, with new shared and source-exclusive components on
every draw. `posterior_epred_y1_and_y2()` returns the marginal means of
the two counts, `mu + lambdaone` and `mu + lambdatwo`. Neither function
conditions on a recorded count. The `y1` slice of
`posterior_predict_y1_and_y2()` has the distribution of the draws of
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
and the `y1` slice of `posterior_epred_y1_and_y2()` equals
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
Both functions return an array in the shape that brms uses for a
multivariate model: posterior draws by observations by response.

The two counts are drawn together: the `y1` and `y2` of each draw share
one draw of the shared component `N_s`, so the covariance of the drawn
counts is `Var(N_s)`, as in the model. The package has no separate
function for `y2`, because draws of `y2` made apart from the draws of
`y1` would not share a draw of `N_s` with them. Marginal draws of `y2`
are `draws[, , "y2"]`, where `draws` is the result of
`posterior_predict_y1_and_y2()`.

## Usage

``` r
posterior_predict_y1_and_y2(
  object,
  newdata = NULL,
  re_formula = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)

posterior_epred_y1_and_y2(
  object,
  newdata = NULL,
  re_formula = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)
```

## Arguments

- object:

  A `brmsfit` from
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
  or
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).

- newdata:

  An optional data frame; if `NULL`, the data of the fit are used. See
  Details.

- re_formula:

  A formula for the group-level effects to include, passed to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html).
  `NULL` (the default) includes all group-level effects, and `NA`
  includes none.

- ndraws, draw_ids:

  Passed to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html)
  to select posterior draws.

- ...:

  Further arguments to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html),
  such as `allow_new_levels`.

## Value

A numeric array with dimensions posterior draws x observations x 2, with
the third dimension named `y1` and `y2`. The values from
`posterior_predict_y1_and_y2()` are integers.

## Details

For each posterior draw of the parameters, the three components `N_s`,
`N1`, and `N2` are drawn independently from their distributions, and
`posterior_predict_y1_and_y2()` returns `y1 = N_s + N1` and
`y2 = N_s + N2`. The marginal moments are

\$\$E\[y_1\] = \mu + \lambda_1, \qquad E\[y_2\] = \mu + \lambda_2,
\qquad \mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_s),\$\$

with \\\mathrm{Var}(N_s) = \mu\\ for `bipois` and \\\mathrm{Var}(N_s) =
\mu + \mu^2 / \phi_s\\ for `binegbin`, where \\\phi_s\\ is `shapes`. The
variance of each count is the sum of the variances of its two
components, for example \\\mathrm{Var}(y_1) = \mu + \mu^2/\phi_s +
\lambda_1 + \lambda_1^2/\phi\_{x1}\\ for `binegbin`.

A predictive check of the pair, such as of the variance of the
difference `y1 - y2`, needs draws of both counts from
`posterior_predict_y1_and_y2()`. Subtracting the recorded `y2` from
draws of
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
treats the two counts as independent and overstates that variance by
\\2\\\mathrm{Var}(N_s)\\.

**New data.** The recorded `y1`, the recorded `y2`, and the observation
indicator do not affect the result. The column `y1` may be omitted from
`newdata`, but brms requires every variable named in the `vint()` term,
so `newdata` must contain `y2` and, for a fit from a `_partialobs()`
constructor, the observation indicator. The values in these columns are
not used, so any integers will do, such as `0` or `-1`. The same is true
for `newdata` in
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).

## See also

[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
for draws of `y1` conditional on the recorded `y2`;
[`vignette("prediction")`](https://anhsmith.github.io/bicountbrms/articles/prediction.md)
for a summary of the prediction functions.

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)

# 500 units. y2 is recorded on every unit, and y1 on about 20% of units.
set.seed(20261115)
dat <- rbicount(
  500,
  mu        = 4,
  lambdaone = 1,
  lambdatwo = 2,
  shapes    = 2,
  shapexone = 4,
  shapextwo = 4
)
dat$y1_obs  <- rbinom(nrow(dat), 1, 0.2)
dat$y1_true <- dat$y1
dat$y1[dat$y1_obs == 0] <- -1L

fit <- brm(
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
  prior    = c(
    prior(normal(1.5, 1), class = "b", nlpar = "logM"),
    prior(normal(0, 1.5), class = "b", nlpar = "logitf"),
    prior(normal(0, 0.5), class = "b", nlpar = "delta"),
    prior(normal(0, 1),   class = "b", nlpar = "kappas",    lb = 0),
    prior(normal(0, 1),   class = "b", nlpar = "kappaxone", lb = 0),
    prior(normal(0, 1),   class = "b", nlpar = "kappaxtwo", lb = 0)
  )
)

# Predictive check of the variance of the difference, on the paired rows
paired <- dat[dat$y1_obs == 1, ]
yrep <- posterior_predict_y1_and_y2(fit, newdata = paired, ndraws = 500)
var_d_rep <- apply(yrep[, , "y1"] - yrep[, , "y2"], 1, var)
hist(var_d_rep)
abline(v = var(paired$y1 - paired$y2))

# Expected counts of the two sources on each paired unit
ep <- posterior_epred_y1_and_y2(fit, newdata = paired)
dim(ep)
} # }
```
