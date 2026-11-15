# Getting started with bicountbrms

``` r

library(brms)
library(bicountbrms)
library(posterior)
library(dplyr)
```

bicountbrms models pairs of counts of the same events recorded by two
sources, such as two observers or two survey methods. The quantities
usually of interest are the amount recorded by the two sources, the
extent to which their records agree, and the direction and size of any
imbalance between them. This vignette simulates paired counts with known
values of these quantities, fits a model that estimates them, and checks
the estimates against the true values.

The package has two brms custom families ([Bürkner
2017](#ref-burkner2017)), `binegbin`, with negative-binomial components,
and `bipois`, with Poisson components, so models are fitted, used for
prediction, and compared with the usual brms functions. Each family has
two constructors:
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
and
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
for data in which both counts were recorded on every row, and
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
for partially paired data, in which source 1 has missing values. This
vignette uses
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md).

## The generative model

Each observed count is the sum of a latent count recorded by both
sources and a latent count recorded by that source alone, a construction
known as ‘trivariate reduction’. Holgate ([1964](#ref-holgate1964)) and
Karlis and Ntzoufras ([2003](#ref-karlis2003)) use this construction
with Poisson components, Kirkpatrick and Neale
([2016](#ref-kirkpatrickApplyingMultivariateDiscrete2016)) and
Kirkpatrick ([2022](#ref-kirkpatrickRMKdiscreteSundryDiscrete2022)) with
negative-binomial components.
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
takes all three latent counts to be independent and negative binomial:

\begin{aligned} N_s &\sim \mathrm{NB2}(\mu,\\ \phi_s) \\ N_1 &\sim
\mathrm{NB2}(\lambda_1,\\ \phi\_{x1}) \qquad N_2 \sim
\mathrm{NB2}(\lambda_2,\\ \phi\_{x2}) \\ y_1 &= N_s + N_1 \qquad y_2 =
N_s + N_2 \end{aligned}

N_s, the shared component, is the latent count recorded by both sources.
N_1 and N_2, the source-exclusive components, are the latent counts
recorded by one source only. The shared component induces the positive
correlation between y_1 and y_2. N_s is never observed, and the
likelihood sums over its possible values k = 0, \dots, \min(y_1, y_2).

`NB2(m, phi)` is `neg_binomial_2` in Stan and
`dnbinom(size = phi, mu = m)` in R: mean `m`, variance `m + m^2/phi`.
Larger `phi` means *less* overdispersion.

In the native parameterisation,
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
has six distributional parameters (dpars), all with a log link:

| Parameter | Meaning | Name in brms |
|----|----|----|
| \mu | rate of the shared component | `mu` |
| \lambda_1, \lambda_2 | rates of the two source-exclusive components | `lambdaone`, `lambdatwo` |
| \phi_s | NB2 dispersion of the shared component | `shapes` |
| \phi\_{x1}, \phi\_{x2} | NB2 dispersions of the two exclusive components | `shapexone`, `shapextwo` |

The dpar `mu` is the rate of the shared component rather than the mean
of either count, so that E\[y_1\] = \mu + \lambda_1. brms requires every
custom family to have a dpar named `mu`.

## The congruence–bias parameterisation

bicountbrms provides an interpretable reparameterisation of the model,
the congruence–bias parameterisation, in which the coordinates directly
represent the quantities usually of interest in a comparison of two
sources: the midpoint M, the average amount recorded by the two sources;
the congruence f, the extent to which their records agree; and the
source bias \beta, the direction and size of any imbalance between them.
The three dispersions are replaced by overdispersions \kappa, for which
\kappa = 0 is a Poisson component. The helper function
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md)
returns the brms formula for this model, and takes one argument for each
coordinate in the table below.

| Coordinate | Meaning | Definition | Name in brms |
|----|----|----|----|
| M | midpoint of the expected counts of the two sources | M = \mu + (\lambda_1 + \lambda_2)/2 | `logM`, for \log M |
| f | congruence, the proportion of M recorded by both sources | f = \mu / M | `logitf`, for \operatorname{logit} f |
| \beta | source bias, zero when the two exclusive rates are equal | \beta = (\lambda_1 - \lambda_2)/(\lambda_1 + \lambda_2) | `delta`, for \operatorname{artanh}\beta |
| \kappa_s, \kappa\_{x1}, \kappa\_{x2} | overdispersions | \kappa = 1/\sqrt{\phi} | `kappas`, `kappaxone`, `kappaxtwo` |

The interactive widget in the article [*The congruence–bias
parameterisation*](https://anhsmith.github.io/bicountbrms/articles/reparameterisation.html)
simulates pairs of counts from one model and shows both
parameterisations side by side: moving a slider for one parameterisation
moves the sliders for the other, and the plot shows the expected
decomposition of each count and the simulated pairs.

We generally recommend the congruence–bias parameterisation. The native
parameterisation might be useful when a predictor is expected to affect
one latent component only, such as a change in the equipment used by one
source (see [Models with predictors](#models-with-predictors)).

## A worked example

### Simulate data from known parameters

The true values are set as midpoint, congruence, source bias, and
overdispersions.
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
converts them to the six native dpars, and
[`rbicount()`](https://anhsmith.github.io/bicountbrms/reference/rbicount.md)
draws the three latent counts and returns them with the two observed
counts.

``` r

set.seed(20260724)

truth <- list(
  M         = 10,
  f         = 0.75,
  beta      = 0.2,
  kappas    = 0.5,
  kappaxone = 0.4,
  kappaxtwo = 0.7
)

dpars <- mfd_to_dpars(
  M         = truth$M,
  f         = truth$f,
  delta     = atanh(truth$beta),
  kappas    = truth$kappas,
  kappaxone = truth$kappaxone,
  kappaxtwo = truth$kappaxtwo
)

dat <- do.call(rbicount, c(list(n = 400), dpars))

head(dat)
#>   y1 y2 Ns N1 N2
#> 1 17 12 11  6  1
#> 2 10 10  7  3  3
#> 3  5  4  2  3  2
#> 4 13 13 10  3  3
#> 5  7 10  5  2  5
#> 6 11 11  9  2  2
cor(dat$y1, dat$y2)
#> [1] 0.8470486
```

The correlation is positive and substantial because both counts include
the same `Ns`.

### A model with intercepts only

A model in the congruence–bias parameterisation is fitted by specifying
its formula with the helper function
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md).
The first argument of
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md)
is the response, `y1 | vint(y2)`: a brms custom family has one response,
so `y1` is the response and `y2` is passed to the family as
supplementary integer data through `vint()`. Each further argument is
the formula of one parameter in the table above, here `~ 1` for all six,
so that each coordinate is the same on every row. The argument
`stanvars = binegbin_stanvars()` adds the Stan function `binegbin_lpmf`
to the `functions` block of the model; without it, the model does not
compile.

``` r

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
  prior    = c(
    prior(normal(2, 1),   class = "b", nlpar = "logM"),
    prior(normal(0, 1.5), class = "b", nlpar = "logitf"),
    prior(normal(0, 0.5), class = "b", nlpar = "delta"),
    prior(normal(0, 1),   class = "b", nlpar = "kappas",    lb = 0),
    prior(normal(0, 1),   class = "b", nlpar = "kappaxone", lb = 0),
    prior(normal(0, 1),   class = "b", nlpar = "kappaxtwo", lb = 0)
  ),
  chains  = 2,
  iter    = 1500,
  refresh = 0,
  backend = "cmdstanr",  # or "rstan"
  seed    = 1
)
```

Unless a prior is specified, brms gives each non-linear parameter a
flat, improper prior, so we recommend setting a prior on every
coordinate. Each prior has `class = "b"` and the name of its parameter.
The prior `normal(2, 1)` on `logM` places 95% of its mass between 1 and
52 on the scale of M, which suits these simulated counts. For other
data, we suggest centring it on the log of a plausible midpoint rather
than widening its SD, because a wider normal prior on the log scale
places more mass on implausibly small and large values ([Smith et al.
2020](#ref-smithInstantaneousVsNoninstantaneous2020), supplement 3, Fig.
S3-1). With an intercept alone in the formula for `logitf`, the prior
`normal(0, 1.5)` gives a density for f within 10% of uniform between 0.1
and 0.9 ([McElreath 2020](#ref-mcelreath2020)). The source bias and the
overdispersions each have a finite null, \delta = 0 for no source bias
and \kappa = 0 for a Poisson component, so their priors have their mode
at zero and shrink the estimates towards the simpler model unless the
data support otherwise ([Simpson et al.
2017](#ref-simpsonPenalisingModelComponent2017)). In the native
parameterisation, a Poisson component is the limit \phi \to \infty,
which no prior on \phi can have as its mode. The bound `lb = 0` makes
each prior on an overdispersion half-normal.

### Convergence

``` r

summary(fit)
#>  Family: binegbin 
#>   Links: mu = log; lambdaone = log; lambdatwo = log; shapes = log; shapexone = log; shapextwo = log 
#> Formula: y1 | vint(y2) ~ 1 
#>          mu ~ logM + log_inv_logit(logitf)
#>          lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)
#>          lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)
#>          shapes ~ -2 * log(kappas)
#>          shapexone ~ -2 * log(kappaxone)
#>          shapextwo ~ -2 * log(kappaxtwo)
#>          logM ~ 1
#>          logitf ~ 1
#>          delta ~ 1
#>          kappas ~ 1
#>          kappaxone ~ 1
#>          kappaxtwo ~ 1
#>    Data: dat (Number of observations: 400) 
#>   Draws: 2 chains, each with iter = 1500; warmup = 750; thin = 1;
#>          total post-warmup draws = 1500
#> 
#> Regression Coefficients:
#>                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> logM_Intercept          2.32      0.03     2.27     2.37 1.00      935      726
#> logitf_Intercept        1.01      0.17     0.67     1.31 1.00      368      560
#> delta_Intercept         0.18      0.03     0.12     0.25 1.00      515      776
#> kappas_Intercept        0.53      0.04     0.46     0.63 1.00      457      729
#> kappaxone_Intercept     0.44      0.10     0.23     0.62 1.00      459      628
#> kappaxtwo_Intercept     0.48      0.17     0.10     0.81 1.00      360      217
#> 
#> Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

`summary(fit)` reports, for each parameter, the potential scale
reduction factor \hat{R} and the bulk and tail effective sample sizes,
and brms warns after sampling if any transitions diverged. A value of
\hat{R} above 1.01 indicates that the chains have not converged to the
same distribution, and an effective sample size below about 100 per
chain indicates too few effectively independent draws for reliable
estimates of the posterior quantiles ([Vehtari et al.
2021](#ref-vehtari2021)). The summary gives each parameter on the scale
of its formula: \log M, \operatorname{logit} f,
\operatorname{artanh}\beta, and \kappa.

### Recovery of the generating values

The table gives the posterior median (`est`) and the 95% credible
interval (`q2.5` to `q97.5`) of each parameter on its natural scale,
with the true value that generated the counts.

``` r

recovery <- as_draws_df(fit) |>
  mutate_variables(
    M         = exp(b_logM_Intercept),
    f         = plogis(b_logitf_Intercept),
    beta      = tanh(b_delta_Intercept),
    kappas    = b_kappas_Intercept,
    kappaxone = b_kappaxone_Intercept,
    kappaxtwo = b_kappaxtwo_Intercept
  ) |>
  subset_draws(variable = names(truth)) |>
  summarise_draws(
    est = median,
    ~ quantile2(.x, probs = c(0.025, 0.975))
  ) |>
  mutate(
    truth   = unlist(truth)[variable],
    covered = truth >= q2.5 & truth <= q97.5
  )

knitr::kable(recovery, digits = 2)
```

| variable  |   est | q2.5 | q97.5 | truth | covered |
|:----------|------:|-----:|------:|------:|:--------|
| M         | 10.22 | 9.71 | 10.72 | 10.00 | TRUE    |
| f         |  0.73 | 0.66 |  0.79 |  0.75 | TRUE    |
| beta      |  0.18 | 0.12 |  0.25 |  0.20 | TRUE    |
| kappas    |  0.53 | 0.46 |  0.63 |  0.50 | TRUE    |
| kappaxone |  0.44 | 0.23 |  0.62 |  0.40 | TRUE    |
| kappaxtwo |  0.48 | 0.10 |  0.81 |  0.70 | TRUE    |

All six 95% intervals contain the true value. The intervals for
`kappaxone` and `kappaxtwo` are the widest, because the overdispersion
of each exclusive component is identified only by the variance of each
count beyond that of the shared component.

### Prediction

[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
predict the response y_1 only, as for any brms family:
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
draws y_1 from its marginal distribution, and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
returns E\[y_1\] = \mu + \lambda_1. y_2 is passed to the model as data
through `vint()`, so neither function returns draws of y_2, and
[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html),
which call these two functions, show y_1 alone.

[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
draws both counts, and returns an array with one row per draw, one
column per row of the data, and one slice for each count. Within each
drawn pair, y_1 and y_2 have the same value of N_s, so the correlation
between the drawn counts is the correlation under the model.
[`posterior_epred_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
returns the expectations of both counts, \mu + \lambda_1 and \mu +
\lambda_2. The table compares the mean, standard deviation, and maximum
of each count, and the correlation between the two counts, with the
averages of the same quantities over 200 drawn datasets.

``` r

yrep <- posterior_predict_y1_and_y2(fit, ndraws = 200)
dim(yrep)
#> [1] 200 400   2

pair_stats <- function(y1, y2) {
  c(
    mean_y1 = mean(y1),
    mean_y2 = mean(y2),
    sd_y1   = sd(y1),
    sd_y2   = sd(y2),
    max_y1  = max(y1),
    max_y2  = max(y2),
    cor     = cor(y1, y2)
  )
}

stats_observed <- pair_stats(dat$y1, dat$y2)

stats_predicted <- yrep |>
  apply(1, \(draw) pair_stats(draw[, "y1"], draw[, "y2"])) |>
  rowMeans()

data.frame(
  observed  = stats_observed,
  predicted = stats_predicted
) |>
  knitr::kable(digits = 2)
```

|         | observed | predicted |
|:--------|---------:|----------:|
| mean_y1 |    10.70 |     10.70 |
| mean_y2 |     9.73 |      9.74 |
| sd_y1   |     5.44 |      5.36 |
| sd_y2   |     5.02 |      5.18 |
| max_y1  |    36.00 |     33.29 |
| max_y2  |    33.00 |     31.89 |
| cor     |     0.85 |      0.84 |

[*Prediction and
imputation*](https://anhsmith.github.io/bicountbrms/articles/prediction.html)
describes the other prediction functions, including
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md),
which draws y_1 conditional on the recorded y_2, and so imputes y_1 on
rows where y_1 was not recorded.

### Model comparison

[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
returns, for each row, the log of the joint probability of (y_1, y_2),
so [`loo()`](https://mc-stan.org/loo/reference/loo.html) estimates the
expected log predictive density (elpd) of the pair ([Vehtari et al.
2017](#ref-vehtari2017)), and models of the same pair can be compared by
their elpd. To compare negative-binomial with Poisson components, we fit
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
which has the same construction as
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
but with Poisson components (and therefore no overdispersions), to the
same data and compare the two fits with
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
Without the three overdispersions,
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md)
returns the formula for
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md).

``` r

fit_pois <- brm(
  bf_mfd(
    y1 | vint(y2),
    logM   ~ 1,
    logitf ~ 1,
    delta  ~ 1
  ),
  family   = bipois(),
  stanvars = bipois_stanvars(),
  data     = dat,
  prior    = c(
    prior(normal(2, 1),   class = "b", nlpar = "logM"),
    prior(normal(0, 1.5), class = "b", nlpar = "logitf"),
    prior(normal(0, 0.5), class = "b", nlpar = "delta")
  ),
  chains  = 2,
  iter    = 1500,
  refresh = 0,
  backend = "cmdstanr",  # or "rstan"
  seed    = 1
)
```

``` r

fit      <- add_criterion(fit, "loo")
fit_pois <- add_criterion(fit_pois, "loo")
loo(fit)
#> 
#> Computed from 1500 by 400 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo  -2176.1 22.5
#> p_loo         6.2  0.6
#> looic      4352.1 45.1
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 1.2]).
#> 
#> All Pareto k estimates are good (k < 0.69).
#> See help('pareto-k-diagnostic') for details.

comparison <- loo_compare(fit, fit_pois)
comparison
#>     model elpd_diff se_diff p_worse diag_diff diag_elpd
#>       fit       0.0     0.0      NA                    
#>  fit_pois    -169.1    30.2    1.00
```

`loo(fit)` prints the elpd estimate of the negative-binomial model and a
table of Pareto k estimates. A row with a Pareto k estimate above about
0.7 has an unreliable estimate of its elpd ([Vehtari et al.
2024](#ref-vehtari2024)).
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
lists the model with the higher elpd first. `elpd_diff` is the
difference in elpd between each model and the first model, and `se_diff`
is the standard error of that difference. The elpd of the Poisson model
is lower by 169 (SE 30), so negative-binomial components describe these
counts much better than Poisson components, as expected from counts
simulated with overdispersions between 0.4 and 0.7.

## Models with predictors

Each coordinate has its own formula. A predictor in the formula for
`logM` multiplies all three rates by the same factor, and a predictor in
the formula for `logitf` or `delta` changes the congruence or the source
bias. The formulas below are for data with a factor `source1_type`, such
as the type of equipment used by source 1, and a grouping factor `group`
with many levels; the simulated data above have neither.

In the model below, the formula for `logM` alone includes a group-level
effect of `group`, so the groups differ in how much the two sources
record but not in the extent to which their records agree. The formula
for `delta` includes `source1_type`, so `source1_type` changes the
balance of exclusive records between the two sources with M and f held
fixed. Each coefficient has a prior with `class = "b"`, as in the model
with intercepts only, and the SD of the group-level effect on `logM` has
a prior with `class = "sd"`.

``` r

fit_groups <- brm(
  bf_mfd(
    y1 | vint(y2),
    logM      ~ 1 + (1 | group),
    logitf    ~ 1,
    delta     ~ source1_type,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  ),
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  data     = my_data,
  prior    = c(
    prior(normal(2, 1),   class = "b",  nlpar = "logM"),
    prior(normal(0, 1),   class = "sd", nlpar = "logM"),
    prior(normal(0, 1.5), class = "b",  nlpar = "logitf"),
    prior(normal(0, 0.5), class = "b",  nlpar = "delta"),
    prior(normal(0, 1),   class = "b",  nlpar = "kappas",    lb = 0),
    prior(normal(0, 1),   class = "b",  nlpar = "kappaxone", lb = 0),
    prior(normal(0, 1),   class = "b",  nlpar = "kappaxtwo", lb = 0)
  ),
  backend = "cmdstanr",  # or "rstan"
  seed    = 1
)
```

``` r

summary(fit_groups)
```

In `summary(fit_groups)`, the row `delta_source1_typeB` is the
difference in \operatorname{artanh}\beta between level `B` of
`source1_type` and the reference level, and the row
`sd(logM_Intercept)`, under `~group`, is the SD of the group-level
effects on \log M.

### Native parameterisation

Where a predictor changes one latent component alone, the native
parameterisation describes it more simply. In the formula below, only
the formula for `lambdaone` includes `source1_type`, so `source1_type`
changes the exclusive rate of source 1 and no other rate, and so changes
M, f, and \beta together. In this parameterisation, the three rates can
be linked across groups only through a group-level effect with a shared
ID in each rate formula, `(1 | s | group)`, under which the three
effects of a group are correlated rather than equal ([Bürkner
2018](#ref-burkner2018)).

``` r

bf(
  y1 | vint(y2) ~ 1 + (1 | s | group),
  lambdaone     ~ 1 + (1 | s | group) + source1_type,
  lambdatwo     ~ 1 + (1 | s | group),
  shapes        ~ 1,
  shapexone     ~ 1,
  shapextwo     ~ 1
)
```

## Further reading

[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
fit partially paired designs, in which y_1 was not recorded on some
rows, with the response `y1 | vint(y2, y1_obs)`. A partially paired
design is fitted in the article [*Partially paired
designs*](https://anhsmith.github.io/bicountbrms/articles/partially-paired-designs.html),
and [*Prediction and
imputation*](https://anhsmith.github.io/bicountbrms/articles/prediction.html)
describes the imputation of the unrecorded values of y_1.

The package `skellambrms` ([Smith
2026](#ref-smithSkellambrmsSkellamDiscreteDifference2026)) models the
difference d = y_1 - y_2 directly. The shared component cancels from d,
so a model of d cannot estimate the midpoint or the congruence.

## References

Bürkner, Paul-Christian. 2017. “brms: An R Package for Bayesian
Multilevel Models Using Stan.” *Journal of Statistical Software* 80 (1):
1–28. <https://doi.org/10.18637/jss.v080.i01>.

Bürkner, Paul-Christian. 2018. “Advanced Bayesian Multilevel Modeling
with the R Package brms.” *The R Journal* 10 (1): 395–411.
<https://doi.org/10.32614/RJ-2018-017>.

Holgate, P. 1964. “Estimation for the Bivariate Poisson Distribution.”
*Biometrika* 51 (1–2): 241–45. <https://doi.org/10.2307/2334210>.

Karlis, Dimitris, and Ioannis Ntzoufras. 2003. “Analysis of Sports Data
by Using Bivariate Poisson Models.” *Journal of the Royal Statistical
Society: Series D (The Statistician)* 52 (3): 381–93.
<https://doi.org/10.1111/1467-9884.00366>.

Kirkpatrick, Robert M. 2022. *RMKdiscrete: Sundry Discrete Probability
Distributions*. <https://CRAN.R-project.org/package=RMKdiscrete>.

Kirkpatrick, Robert M., and Michael C. Neale. 2016. “Applying
Multivariate Discrete Distributions to Genetically Informative Count
Data.” *Behavior Genetics* 46 (2): 252–68.
<https://doi.org/10.1007/s10519-015-9757-z>.

McElreath, Richard. 2020. *Statistical Rethinking: A Bayesian Course
with Examples in R and Stan*. 2nd ed. Chapman; Hall/CRC.
<https://doi.org/10.1201/9780429029608>.

Simpson, Daniel, Håvard Rue, Andrea Riebler, Thiago G. Martins, and
Sigrunn H. Sørbye. 2017. “Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors.” *Statistical
Science* 32 (1): 1–28. <https://doi.org/10.1214/16-STS576>.

Smith, A. N. H., D. Acuña-Marrero, P. Salinas-de-León, E. S. Harvey, M.
D. M. Pawley, and M. J. Anderson. 2020. “Instantaneous Vs.
Non-Instantaneous Diver-Operated Stereo-Video (DOV) Surveys of Highly
Mobile Sharks in the Galápagos Marine Reserve.” *Marine Ecology Progress
Series* 649: 111–23. <https://doi.org/10.3354/meps13447>.

Smith, Adam N. H. 2026. *skellambrms: Skellam and Discrete-Difference
brms Custom Families for Count Differences*.
<https://doi.org/10.5281/zenodo.22231870>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved \widehat{R} for Assessing Convergence of MCMC
(with Discussion).” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2024. “Pareto Smoothed Importance Sampling.” *Journal of Machine
Learning Research* 25 (72): 1–58.
