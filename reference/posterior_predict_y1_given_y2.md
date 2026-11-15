# Draws and expectations of y1 given y2

`posterior_predict_y1_given_y2()` draws `y1` from its distribution
conditional on the recorded `y2` and the model parameters, \\p(y_1 \mid
y_2, \theta)\\. `posterior_epred_y1_given_y2()` returns the conditional
expectation \\E\[y_1 \mid y_2\]\\. The main use of both functions is to
impute `y1` on the rows where `y1` was not recorded. On an unpaired row
of a fit from
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
or
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md),
the draws are the imputed values of the unrecorded `y1`, and the
expectation is the mean of those imputed values.

[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
do not condition on `y2`. For a fit from any of the four constructors,
these two functions return draws from the marginal distribution of `y1`
and its marginal mean `mu + lambdaone`, as for any other brms family.
[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
draws both counts.

## Usage

``` r
posterior_predict_y1_given_y2(
  object,
  newdata = NULL,
  re_formula = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  keep_observed = FALSE,
  ...
)

posterior_epred_y1_given_y2(
  object,
  newdata = NULL,
  re_formula = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  keep_observed = FALSE,
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

- keep_observed:

  Logical. If `FALSE` (the default), every row returns a conditional
  draw or the conditional expectation. If `TRUE`, paired rows return the
  recorded `y1`.

- ...:

  Further arguments to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html),
  such as `allow_new_levels`.

## Value

A numeric matrix with posterior draws in rows and observations in
columns. The values from `posterior_predict_y1_given_y2()` are integers.

## Details

Both families construct the pair as `y1 = N_s + N1` and `y2 = N_s + N2`,
where the shared component `N_s` has mean `mu` and dispersion `shapes`,
and the source-exclusive components `N1` and `N2` have means `lambdaone`
and `lambdatwo` and dispersions `shapexone` and `shapextwo`. The three
components are independent given the parameters, so `N1` is independent
of `y2`, and \\y_1 \mid y_2\\ is the sum of \\N_s \mid y_2\\ and an
independent draw of `N1`. The conditional distribution of the shared
component is

\$\$P(N_s = k \mid y_2) \propto p(k \mid \mu, \phi_s)\\ p(y_2 - k \mid
\lambda_2, \phi\_{x2}), \quad k = 0, \ldots, y_2,\$\$

where \\p(\cdot \mid m, \phi)\\ is the negative-binomial probability
mass function with mean \\m\\ and inverse dispersion \\\phi\\ for
`binegbin`, and the Poisson probability mass function with mean \\m\\
for `bipois`. For `bipois`, this distribution is Binomial with size
\\y_2\\ and probability \\\mu / (\mu + \lambda_2)\\. For `binegbin`, the
distribution has no closed form and is evaluated over the \\y_2 + 1\\
values of \\k\\. In both families,

\$\$E\[y_1 \mid y_2\] = E\[N_s \mid y_2\] + \lambda_1.\$\$

The exclusive-component dispersion in the conditional is `shapextwo`,
the dispersion of the second exclusive component, because the count
conditioned on is `y2`.

**Recorded values of `y1`.** With `keep_observed = FALSE`, the default,
every row receives a conditional draw or the conditional expectation,
including the paired rows. On a paired row, the conditional draw is the
quantity to compare with the recorded `y1` in a conditional predictive
check. With `keep_observed = TRUE`, every draw on a paired row is the
recorded `y1`, and the draws on unpaired rows are unchanged. Summed over
rows, the draws returned with `keep_observed = TRUE` are draws of the
total of `y1`: recorded where `y1` was recorded, and imputed where `y1`
was not recorded. To reproduce results from
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
in bicountbrms 0.10.1, use these two functions with
`keep_observed = FALSE`; given the same seed, the draws are identical.

**New data.** `newdata` must contain the recorded `y2` of each row,
because the result depends on `y2`. For a fit from a `_partialobs()`
constructor, brms also requires the observation indicator named in
`vint()`. The column `y1` may be omitted, except that with
`keep_observed = TRUE` a row indicated as paired needs its recorded
`y1`, and the function stops with an error if that `y1` is missing.

Stored fits whose family object declares the pre-0.7.0 rate names
(`lambdaem`, `lambdalb`) or a single exclusive-component dispersion
(`shapex`) are read under those names, as by the brms methods.

## See also

[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
for draws of both counts from the marginal model;
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
  bf(
    y1 | vint(y2, y1_obs) ~ 1,
    lambdaone ~ 1,
    lambdatwo ~ 1,
    shapes    ~ 1,
    shapexone ~ 1,
    shapextwo ~ 1
  ),
  family   = binegbin_partialobs(),
  stanvars = binegbin_partialobs_stanvars(),
  data     = dat
)

# Draws of y1 on every unit: recorded where y1_obs is 1,
# imputed from the y2 of each unit where y1_obs is 0.
y1_draws <- posterior_predict_y1_given_y2(fit, keep_observed = TRUE)

# Posterior of the source-1 total, against the simulated total
quantile(rowSums(y1_draws), probs = c(0.025, 0.5, 0.975))
sum(dat$y1_true)

# Imputed draws and E[y1 | y2] on the units where source 1 was not
# recorded only
unrec <- dat$y1_obs == 0
y1_imp <- posterior_predict_y1_given_y2(fit)[, unrec]
ep_imp <- posterior_epred_y1_given_y2(fit)[, unrec]
} # }
```
