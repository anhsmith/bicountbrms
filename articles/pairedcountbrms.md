# Getting started with pairedcountbrms

``` r

library(brms)
#> Warning: package 'Rcpp' was built under R version 4.6.1
library(pairedcountbrms)
```

This vignette walks the package’s R API end to end on one family,
[`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md):
simulate a matched pair of counts from known parameters, fit them, check
that the posterior finds the truth, and then predict and score. The
other joint family
([`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md))
works identically with two fewer dpars; the difference families
([`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md),
[`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md),
and so on) are shown in the README.

Everything here is a brms custom family ([Bürkner
2017](#ref-burkner2017)), so the fitting, prediction and
model-comparison interfaces are brms’s own. The non-linear formula
syntax used further down is documented in Bürkner
([2018](#ref-burkner2018)).

## The generative model

[`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
builds a bivariate count pair by **trivariate reduction** ([Holgate
1964](#ref-holgate1964); [Karlis and Ntzoufras 2003](#ref-karlis2003)).
Three independent Negative-Binomial counts are drawn, and the two
observed counts share one of them:

``` math
\begin{aligned}
N_{\text{shared}} &\sim \mathrm{NB2}(\mu,\ \phi_s) \\
N_1 &\sim \mathrm{NB2}(\lambda_1,\ \phi_x) \qquad
N_2 \sim \mathrm{NB2}(\lambda_2,\ \phi_x) \\[4pt]
y_1 &= N_{\text{shared}} + N_1 \qquad
y_2 = N_{\text{shared}} + N_2
\end{aligned}
```

The shared component is what both sources saw; the two excesses are what
each saw alone. $`N_{\text{shared}}`$ is never observed — it is
marginalised out analytically in the likelihood — but it is what induces
the correlation between the pair.

`NB2(m, phi)` is Stan’s `neg_binomial_2` and R’s
`dnbinom(size = phi, mu = m)`: mean `m`, variance `m + m^2/phi`. Larger
`phi` means *less* overdispersion.

Five dpars, all with a log link:

| dpar        | role                                              |
|-------------|---------------------------------------------------|
| `mu`        | rate of the shared component, $`\mu`$             |
| `lambdaone` | rate of the first source’s excess, $`\lambda_1`$  |
| `lambdatwo` | rate of the second source’s excess, $`\lambda_2`$ |
| `shapes`    | dispersion $`\phi_s`$ of the shared component     |
| `shapex`    | dispersion $`\phi_x`$, shared by both excesses    |

`mu` is brms’s mandatory first-dpar name. Here it is bound to the shared
component’s *rate* — it is not the mean of either response, and not the
mean of their difference. `E[y1] = mu + lambdaone`.

The two rates are spelled `lambdaone`/`lambdatwo` rather than
`lambda1`/`lambda2` because
[`custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
rejects dpar names ending in a digit. The maths below writes them
$`\lambda_1`$ and $`\lambda_2`$.

All five use a log link, the conventional log-linear rate
parameterisation for this construction ([Karlis and Ntzoufras
2003](#ref-karlis2003)).

## Simulate from known parameters

``` r

set.seed(20260724)

n <- 400
truth <- list(
  mu        = 8,
  lambdaone = 3,
  lambdatwo = 2,
  shapes    = 4,
  shapex    = 6
)

n_shared <- rnbinom(n, size = truth$shapes, mu = truth$mu)
n1       <- rnbinom(n, size = truth$shapex, mu = truth$lambdaone)
n2       <- rnbinom(n, size = truth$shapex, mu = truth$lambdatwo)

dat <- data.frame(
  y1 = n_shared + n1,
  y2 = n_shared + n2
)

str(dat)
#> 'data.frame':    400 obs. of  2 variables:
#>  $ y1: num  11 8 7 15 5 13 9 15 13 8 ...
#>  $ y2: num  11 7 6 16 6 11 9 15 9 7 ...
cor(dat$y1, dat$y2)
#> [1] 0.8911927
```

The correlation is positive and substantial because both counts carry
the same `n_shared`. That shared term is exactly what a model of the
difference alone throws away.

## Fit

Two things are specific to this package and easy to miss.

**The second count travels via `vint()`.**
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
declares a single response column, so only `y1` can be the response.
`y2` is passed alongside as supplementary integer data, and the family’s
Stan signature reads it from there.

**`stanvars` is not optional.**
[`binegbin_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
injects the `binegbin_lpmf` Stan function into the model’s `functions`
block. Without it the generated model will not compile.

``` r

fit <- brm(
  bf(
    y1 | vint(y2) ~ 1,
    lambdaone     ~ 1,
    lambdatwo     ~ 1,
    shapes        ~ 1,
    shapex        ~ 1
  ),
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  data     = dat,
  prior    = c(
    prior(normal(2, 1), class = "Intercept"),
    prior(normal(2, 1), class = "Intercept", dpar = "lambdaone"),
    prior(normal(2, 1), class = "Intercept", dpar = "lambdatwo"),
    prior(normal(2, 1), class = "Intercept", dpar = "shapes"),
    prior(normal(2, 1), class = "Intercept", dpar = "shapex")
  ),
  chains  = 2,
  iter    = 1500,
  refresh = 0,
  backend = backend,
  seed    = 1
)
```

One weakly-informative prior serves all five, since every dpar here is a
log-linked positive quantity of similar magnitude. On the natural scale
`normal(2, 1)` is lognormal with 95% of its mass in roughly `[1, 52]`.

The five truths span `0.69` to `2.08` on the log scale, so no single
prior sits away from all of them: the shared rate’s truth falls near the
prior mean, the second excess rate’s 1.3 prior SDs below it. Recovery
here is therefore not on its own evidence that the prior is
uninfluential. The article [*The anatomy of a paired
count*](https://anhsmith.github.io/pairedcountbrms/articles/paired-count-anatomy.html)
tests that directly, shifting a prior median sevenfold and moving the
corresponding posterior median by 0.005.

To adapt these priors, shift the **mean** to the scale of the counts
rather than increasing the SD. All five dpars are log-linked, so a
normal prior is lognormal on the natural scale, and increasing its SD
moves mass to implausible values rather than making the prior neutral
([Smith et al. 2020](#ref-smithInstantaneousVsNoninstantaneous2020),
supplement 3).

## Did it recover the truth?

All five dpars are log-linked, so each posterior intercept exponentiates
back onto the natural scale.

``` r

draws <- as_draws_df(fit)

pars <- c(
  mu        = "b_Intercept",
  lambdaone = "b_lambdaone_Intercept",
  lambdatwo = "b_lambdatwo_Intercept",
  shapes    = "b_shapes_Intercept",
  shapex    = "b_shapex_Intercept"
)

recovery <- do.call(rbind, lapply(names(pars), function(p) {
  x <- exp(draws[[pars[[p]]]])
  data.frame(
    dpar  = p,
    truth = truth[[p]],
    est   = median(x),
    lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975))
  )
}))
recovery$covered <- with(recovery, truth >= lower & truth <= upper)

knitr::kable(recovery, digits = 2)
```

| dpar      | truth |  est | lower | upper | covered |
|:----------|------:|-----:|------:|------:|:--------|
| mu        |     8 | 8.36 |  7.65 |  9.04 | TRUE    |
| lambdaone |     3 | 2.87 |  2.39 |  3.41 | TRUE    |
| lambdatwo |     2 | 1.95 |  1.47 |  2.49 | TRUE    |
| shapes    |     4 | 3.82 |  2.97 |  4.95 | TRUE    |
| shapex    |     6 | 7.10 |  3.10 | 23.27 | TRUE    |

The three rates should land tightly on their true values. The two
dispersions are estimated from an aggregate mean–variance mismatch
rather than from any directly observed quantity, so their intervals are
wider — `shapex` especially, since it is identified only through the
part of the pair’s spread that the shared component cannot explain. Wide
but covering is the expected result here, not a warning sign.

## Predict and score

[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
draws new `y1` **conditional on each row’s observed `y2`**. It samples
the discrete conditional distribution of $`N_{\text{shared}}`$ given
`y2`, then adds a fresh excess draw — the exact conditional, not an
approximation.

``` r

yrep <- posterior_predict(fit, ndraws = 200)
dim(yrep)
#> [1] 200 400

data.frame(
  quantity = c("mean", "sd", "max"),
  observed = c(mean(dat$y1), sd(dat$y1), max(dat$y1)),
  predicted = c(mean(yrep), mean(apply(yrep, 1, sd)), mean(apply(yrep, 1, max)))
) |>
  knitr::kable(digits = 2)
```

| quantity | observed | predicted |
|:---------|---------:|----------:|
| mean     |    11.25 |     11.24 |
| sd       |     5.59 |      5.49 |
| max      |    37.00 |     34.12 |

[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
gives the pointwise log-likelihood of the *joint* pair, evaluated by an
independent R implementation of the same marginalisation sum that the
Stan function computes. It feeds
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`waic()`](https://mc-stan.org/loo/reference/waic.html) in the usual way
([Vehtari et al. 2017](#ref-vehtari2017)). The Pareto $`k`$ diagnostic
reported alongside the estimate flags observations whose importance
weights are unreliable ([Vehtari et al. 2024](#ref-vehtari2024)).

``` r

ll <- log_lik(fit)
dim(ll)
#> [1] 1500  400

loo(fit)
#> 
#> Computed from 1500 by 400 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo  -2146.3 22.7
#> p_loo         4.8  0.5
#> looic      4292.7 45.4
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.3]).
#> 
#> All Pareto k estimates are good (k < 0.69).
#> See help('pareto-k-diagnostic') for details.
```

## One limitation to know about

[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
— and its aliases
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
— **errors on any truncated custom-family fit** in brms 2.23.0. The
dispatcher checks for truncation before it checks the family type, and
the truncated branch has no fallback to `posterior_epred_custom()`. This
is a brms limitation, not a bug in this package.

On a truncated fit, call the family’s method directly:

``` r

posterior_epred_binegbin(prepare_predictions(fit))
```

[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
are unaffected and work correctly with truncation.

## Where to go next

- [`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  — the non-overdispersed Poisson sibling. Same construction, three
  dpars instead of five. Use it when the margins are not overdispersed;
  compare the two with
  [`loo()`](https://mc-stan.org/loo/reference/loo.html).
- [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
  — censoring-aware, for data where one of the two counts is missing on
  some rows. It admits those rows via the integrated-out marginal
  instead of dropping them. It also gives each excess component its own
  dispersion (`shapexone`/`shapextwo`, six dpars rather than five), so
  the two sources may be differently overdispersed; tying them together
  with
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  recovers the single-`shapex` model shown above.
- The $`(M, f, \delta)`$**reparameterisation** — overall level $`M`$,
  congruence $`f`$, and source bias $`\delta`$, with the dispersions on
  an SD scale $`\kappa = 1/\sqrt{\phi}`$ where $`\kappa = 0`$ is the
  Poisson limit. It is fitted through
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  rather than a separate family, and
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  converts between the two. Because $`\delta = 0`$ and $`\kappa = 0`$
  are finite, interpretable nulls, these coordinates admit
  penalised-complexity priors ([Simpson et al.
  2017](#ref-simpsonPenalisingModelComponent2017)), which shrink to the
  simpler model unless the data support otherwise; the rate
  parameterisation has no such nulls. Worked end to end in the article
  [*The anatomy of a paired
  count*](https://anhsmith.github.io/pairedcountbrms/articles/paired-count-anatomy.html).
- [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  /
  [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md),
  [`dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  /
  [`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md),
  [`dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  /
  [`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  — the difference families, which model `d = y1 - y2` directly. The `1`
  variants fix the location at zero (does the pair agree on average?);
  the `2` variants estimate it (by how much do they disagree?). All
  support truncation through
  [`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html).
  The difference of two independent Poisson counts is
  Skellam-distributed ([Skellam 1946](#ref-skellam1946)); for the
  Bayesian treatment of count differences generally, see Karlis and
  Ntzoufras ([2006](#ref-karlis2006)).

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

Karlis, Dimitris, and Ioannis Ntzoufras. 2006. “Bayesian Analysis of the
Differences of Count Data.” *Statistics in Medicine* 25 (11): 1885–905.
<https://doi.org/10.1002/sim.2382>.

Simpson, Daniel, Håvard Rue, Andrea Riebler, Thiago G. Martins, and
Sigrunn H. Sørbye. 2017. “Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors.” *Statistical
Science* 32 (1): 1–28. <https://doi.org/10.1214/16-STS576>.

Skellam, J. G. 1946. “The Frequency Distribution of the Difference
Between Two Poisson Variates Belonging to Different Populations.”
*Journal of the Royal Statistical Society* 109 (3): 296.
<https://doi.org/10.2307/2981372>.

Smith, A. N. H., D. Acuña-Marrero, P. Salinas-de-León, E. S. Harvey, M.
D. M. Pawley, and M. J. Anderson. 2020. “Instantaneous Vs.
Non-Instantaneous Diver-Operated Stereo-Video (DOV) Surveys of Highly
Mobile Sharks in the Galápagos Marine Reserve.” *Marine Ecology Progress
Series* 649: 111–23. <https://doi.org/10.3354/meps13447>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

Vehtari, Aki, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
Gabry. 2024. “Pareto Smoothed Importance Sampling.” *Journal of Machine
Learning Research* 25 (72): 1–58.
