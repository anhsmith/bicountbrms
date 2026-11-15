# Prediction and imputation

Both families model a pair of counts as y_1 = N_s + N_1 and y_2 = N_s +
N_2, where N_s is the shared component, recorded by both sources, and
N_1 and N_2 are the source-exclusive components. The three components
are independent given their means \mu (`mu`), \lambda_1 (`lambdaone`),
and \lambda_2 (`lambdatwo`). Six functions return predictions from a
fit, either as draws (the `posterior_predict*()` functions) or as
expectations (the `posterior_epred*()` functions). In the table, \theta
denotes the parameters under one posterior draw.

| Function | Distribution or quantity | Returns |
|----|----|----|
| [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html) | p(y_1 \mid \theta) | draws × observations |
| [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html) | E\[y_1\] = \mu + \lambda_1 | draws × observations |
| [`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md) | p(y_1, y_2 \mid \theta) | draws × observations × 2 |
| [`posterior_epred_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md) | E\[y_1\] = \mu + \lambda_1, E\[y_2\] = \mu + \lambda_2 | draws × observations × 2 |
| [`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md) | p(y_1 \mid y_2, \theta) | draws × observations |
| [`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md) | E\[y_1 \mid y_2\] = E\[N_s \mid y_2\] + \lambda_1 | draws × observations |

The third dimension of the array returned by the two `_y1_and_y2()`
functions is named `y1` and `y2`:

``` r

draws <- posterior_predict_y1_and_y2(fit)
draws[, , "y1"]  # one row per posterior draw, one column per observation
draws[, , "y2"]
```

The first four functions in the table depend on the parameters alone.
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
and
[`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
also condition on the recorded y_2, and are used mainly to impute y_1 on
the unpaired rows of a fit from
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
or
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md),
where y_1 was not recorded. With `keep_observed = TRUE`,
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
and
[`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
return the recorded y_1 on paired rows.

## brms functions for post-processing

| Function | Quantity |
|----|----|
| [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`predict()`](https://rdrr.io/r/stats/predict.html) | Marginal mean and draws of y_1, summarised |
| [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html) | Marginal mean of y_1 over a grid of predictor values |
| [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html) | Recorded y_1 against marginal draws of y_1 |
| [`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html), [`loo_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html) | Variance of the marginal mean of y_1, divided by that variance plus the residual variance |
| LOO-PIT, through `pp_check(type = "loo_pit_overlay")` | Leave-one-out probability integral transform of the recorded y_1 |
| [`loo()`](https://mc-stan.org/loo/reference/loo.html), [`waic()`](https://mc-stan.org/loo/reference/waic.html), and [`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html) | Joint probability of both counts on paired rows, and probability of y_2 on unpaired rows |

On an unpaired row of a `_partialobs()` fit, `y1` is not a recorded
count.
[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
[`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html),
[`loo_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html),
and LOO-PIT compare the predictions with y_1. For such a fit, these four
functions are meaningful only with `newdata` restricted to the paired
rows, as in `bayes_R2(fit, newdata = subset(dat, y1_obs == 1))`.
[`loo()`](https://mc-stan.org/loo/reference/loo.html),
[`waic()`](https://mc-stan.org/loo/reference/waic.html), and
[`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html) need
no such restriction, because the likelihood term for an unpaired row
does not depend on y_1.

## Conditional and joint draws

Given y_2, the shared component has the distribution

P(N_s = k \mid y_2) \propto p(k \mid \mu, \phi_s)\\ p(y_2 - k \mid
\lambda_2, \phi\_{x2}), \qquad k = 0, \ldots, y_2,

where p(\cdot \mid m, \phi) is the negative-binomial probability mass
function with mean m and (inverse) dispersion parameter \phi for
`binegbin`, and the Poisson probability mass function with mean m for
`bipois`. For `bipois`, this distribution is Binomial with size y_2 and
probability \mu / (\mu + \lambda_2).
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
adds an independent draw of N_1 to a draw of N_s from this distribution.

[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
draws N_s once for each pair and adds the same draw to both counts, so
the drawn counts have \mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_s), as in
the model. The package has no separate function for y_2, because draws
of y_2 made separately would not share their draws of N_s with the draws
of y_1.

## Imputing missing values of y_1

To impute missing values of y_1 conditional on the recorded values of
y_2:

1.  Code the observation indicator `y1_obs` as `1` on paired rows and
    `0` on unpaired rows, and set `y1` to any integer on the unpaired
    rows. We suggest `-1`, because it cannot be confused with a recorded
    count. `y1` must not be `NA`.
2.  Fit the model with
    [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
    (or
    [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)),
    passing the indicator as the second argument of `vint()`.
3.  Draw the unrecorded counts with
    [`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md),
    which draws each y_1 given the y_2 recorded on its row.

``` r

# 500 units. y2 is recorded on every unit, and y1 on each unit with
# probability 0.2.
set.seed(20261115)
dpars <- mfd_to_dpars(
  M      = 5,
  f      = 0.7,
  delta  = atanh(0.2),
  kappas = 0.7,
  kappax = 0.5
)
dat <- do.call(rbicount, c(list(n = 500), dpars))
dat$y1_true <- dat$y1
dat$y1[rbinom(nrow(dat), 1, 0.2) == 0] <- NA

# Add the indicator y1_obs, 1 where y1 was recorded, and set the missing values
# of y1 to -1, so that brms keeps those rows
dat$y1_obs <- as.integer(!is.na(dat$y1))
dat$y1[is.na(dat$y1)] <- -1L

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

# Draws of the unrecorded counts: one column per unpaired row
y1_imputed <- posterior_predict_y1_given_y2(fit)[, dat$y1_obs == 0]

# Draws of y1 on every row: recorded values on paired rows, imputed values
# on unpaired rows
y1_completed <- posterior_predict_y1_given_y2(fit, keep_observed = TRUE)
# Posterior of the total of y1 over all rows, against the simulated total
quantile(rowSums(y1_completed), c(0.025, 0.5, 0.975))
sum(dat$y1_true)
```

Each row of `y1_imputed` is one posterior draw of the unrecorded counts,
so summaries of the draws include both the uncertainty in the parameters
and the variation of the counts around their expected values.
[`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
gives the conditional expectation E\[y_1 \mid y_2\] instead, which omits
the variation of the counts. In the article [*Partially paired
designs*](https://anhsmith.github.io/bicountbrms/articles/partially-paired-designs.html),
imputed values of y_1 are compared with values withheld from the model.

## Checking the model against both counts

``` r

# Variance of d = y1 - y2 on the paired rows, against marginal draws
paired <- subset(dat, y1_obs == 1)
yrep   <- posterior_predict_y1_and_y2(fit, newdata = paired, ndraws = 500)
var_d  <- apply(yrep[, , "y1"] - yrep[, , "y2"], 1, var)
mean(var_d >= var(paired$y1 - paired$y2))
```

The last line is the posterior predictive p-value for the variance of d.
A value near 0 or 1 suggests that the model misstates the spread of the
difference.
