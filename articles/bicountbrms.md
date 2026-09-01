# Getting started with bicountbrms

``` r

library(brms)
library(bicountbrms)
```

This vignette demonstrates how to use one family,
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
from simulation to prediction: counts are simulated from known
parameters, fitted, and the credible intervals checked against the
values that generated them. The other joint families
([`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
and the partially observed
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md))
follow the same pattern; the README compares them.

Everything here is a brms custom family ([Bürkner
2017](#ref-burkner2017)), so the fitting, prediction and
model-comparison interfaces are brms’s own. The non-linear formula
syntax used further down is documented in Bürkner
([2018](#ref-burkner2018)).

## The generative model

Each observed count is the sum of a latent count that both sources
recorded and one that only that source recorded. This construction is
**trivariate reduction**: Holgate ([1964](#ref-holgate1964)) and Karlis
and Ntzoufras ([2003](#ref-karlis2003)) use it with Poisson components,
Kirkpatrick and Neale
([2016](#ref-kirkpatrickApplyingMultivariateDiscrete2016)) and
Kirkpatrick ([2022](#ref-kirkpatrickRMKdiscreteSundryDiscrete2022)) with
negative-binomial components.
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
takes all three latent counts to be independent and negative binomial:

``` math
\begin{aligned}
N_{\text{shared}} &\sim \mathrm{NB2}(\mu,\ \phi_s) \\
N_1 &\sim \mathrm{NB2}(\lambda_1,\ \phi_{x1}) \qquad
N_2 \sim \mathrm{NB2}(\lambda_2,\ \phi_{x2}) \\[4pt]
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

Six dpars, all with a log link:

| dpar        | role                                                   |
|-------------|--------------------------------------------------------|
| `mu`        | rate of the shared component, $`\mu`$                  |
| `lambdaone` | rate of the first source’s excess, $`\lambda_1`$       |
| `lambdatwo` | rate of the second source’s excess, $`\lambda_2`$      |
| `shapes`    | dispersion $`\phi_s`$ of the shared component          |
| `shapexone` | dispersion $`\phi_{x1}`$ of the first source’s excess  |
| `shapextwo` | dispersion $`\phi_{x2}`$ of the second source’s excess |

[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
rejects any family whose dpars do not include one named literally `mu`.
The requirement is on the name, not on the position: `mu` may appear
anywhere in the vector. Here it is bound to the shared component’s
*rate*, so it is not the mean of either response, nor the mean of their
difference. `E[y1] = mu + lambdaone`.

The two rates are spelled `lambdaone`/`lambdatwo` rather than
`lambda1`/`lambda2` because
[`custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
rejects dpar names ending in a digit. The maths below writes them
$`\lambda_1`$ and $`\lambda_2`$.

All six use a log link, the conventional log-linear rate
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
  shapexone = 6,
  shapextwo = 2
)

n_shared <- rnbinom(n, size = truth$shapes, mu = truth$mu)
n1       <- rnbinom(n, size = truth$shapexone, mu = truth$lambdaone)
n2       <- rnbinom(n, size = truth$shapextwo, mu = truth$lambdatwo)

dat <- data.frame(
  y1 = n_shared + n1,
  y2 = n_shared + n2
)

str(dat)
#> 'data.frame':    400 obs. of  2 variables:
#>  $ y1: num  11 8 7 15 5 13 9 15 13 8 ...
#>  $ y2: num  11 7 5 16 6 11 9 16 9 7 ...
cor(dat$y1, dat$y2)
#> [1] 0.8518386
```

The correlation is positive and substantial because both counts include
the same `n_shared`.

## Fit

Two things are specific to this package.

**The second count is supplied through `vint()`.**
[`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
declares a single response column, so only `y1` can be the response.
`y2` is passed alongside as supplementary integer data, and the family’s
Stan signature reads it from there.

**`stanvars` is not optional.**
[`binegbin_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
injects the `binegbin_lpmf` Stan function into the model’s `functions`
block. Without it, the generated model will not compile.

``` r

fit <- brm(
  bf(
    y1 | vint(y2) ~ 1,
    lambdaone     ~ 1,
    lambdatwo     ~ 1,
    shapes        ~ 1,
    shapexone     ~ 1,
    shapextwo     ~ 1
  ),
  family   = binegbin(),
  stanvars = binegbin_stanvars(),
  data     = dat,
  prior    = c(
    prior(normal(2, 1), class = "Intercept"),
    prior(normal(2, 1), class = "Intercept", dpar = "lambdaone"),
    prior(normal(2, 1), class = "Intercept", dpar = "lambdatwo"),
    prior(normal(2, 1), class = "Intercept", dpar = "shapes"),
    prior(normal(2, 1), class = "Intercept", dpar = "shapexone"),
    prior(normal(2, 1), class = "Intercept", dpar = "shapextwo")
  ),
  chains  = 2,
  iter    = 1500,
  refresh = 0,
  backend = backend,
  seed    = 1
)
```

One weakly-informative prior serves all six, since **every dpar here is
a log-linked positive quantity** of similar magnitude. On the natural
scale `normal(2, 1)` is lognormal, with 95% of its mass between 1 and
52.

The six truths span 0.69 to 2.08 on the log scale, so no single prior
sits away from all of them: the shared rate’s truth falls near the prior
mean, the second excess rate’s 1.3 prior SDs below it. Recovery here is
therefore not on its own evidence that the prior is uninfluential. The
article [*The anatomy of a paired
count*](https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.html)
tests that directly, shifting a prior median sevenfold and moving the
corresponding posterior median by 0.005.

To adapt these priors, shift the **mean** to the scale of the counts
rather than increasing the SD. All six dpars are log-linked, so a normal
prior is lognormal on the natural scale, and increasing its SD moves
mass to implausible values rather than making the prior neutral ([Smith
et al. 2020](#ref-smithInstantaneousVsNoninstantaneous2020), supplement
3, Fig. S3-1).

## Recovery of the generating values

All six dpars are log-linked, so each posterior intercept exponentiates
back onto the natural scale.

``` r

draws <- as_draws_df(fit)

pars <- c(
  mu        = "b_Intercept",
  lambdaone = "b_lambdaone_Intercept",
  lambdatwo = "b_lambdatwo_Intercept",
  shapes    = "b_shapes_Intercept",
  shapexone = "b_shapexone_Intercept",
  shapextwo = "b_shapextwo_Intercept"
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
| mu        |     8 | 7.91 |  7.14 |  8.61 | TRUE    |
| lambdaone |     3 | 3.33 |  2.81 |  3.97 | TRUE    |
| lambdatwo |     2 | 2.29 |  1.77 |  2.94 | TRUE    |
| shapes    |     4 | 3.58 |  2.68 |  4.73 | TRUE    |
| shapexone |     6 | 7.16 |  3.31 | 21.90 | TRUE    |
| shapextwo |     2 | 3.13 |  1.44 | 10.05 | TRUE    |

The three rates should land tightly on their true values. The three
dispersions are estimated from an aggregate mean–variance mismatch
rather than from any directly observed quantity, so their intervals are
wider; the two excess dispersions especially, since they are identified
only through the part of the pair’s spread that the shared component
cannot explain. Wide but covering is the expected result here, not a
warning sign.

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
| sd       |     5.59 |      5.43 |
| max      |    37.00 |     33.23 |

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
#> elpd_loo  -2202.8 22.6
#> p_loo         5.8  0.6
#> looic      4405.7 45.2
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 1.5]).
#> 
#> All Pareto k estimates are good (k < 0.69).
#> See help('pareto-k-diagnostic') for details.
```

## Quantities returned by `posterior_epred()`

Both families return the same quantity:

``` math
\mathrm{E}[y_1 \mid y_2] = \mathrm{E}[N_{\text{shared}} \mid y_2] + \lambda_1,
```

the expected first count *given the second one that was actually
observed*, not a marginal expectation. It is what the second source
implies about the first.

For
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
$`\mathrm{E}[N_{\text{shared}} \mid y_2]`$ is closed form, because
conditioning a sum of independent Poissons on its total gives a Binomial
split. For
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
there is no such shortcut, so the conditional is evaluated over its
support $`k = 0, \ldots, y_2`$ and summed. Both are exact, and each
agrees with the mean of its own
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
draws to within Monte Carlo error.

The conditional expectation is returned for every row, including — under
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
— those where $`y_1`$ was never observed. Imputing the unobserved margin
is what those two constructors are for.

`posterior_epred(fit)` needs no special handling: it dispatches to the
family’s method, as do
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html).

## Related families and further reading

- [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  — the non-overdispersed Poisson sibling. Same construction, three
  dpars instead of six. Use it when the margins are not overdispersed;
  compare the two with
  [`loo()`](https://mc-stan.org/loo/reference/loo.html).
- **The symmetric model.** The fit above let the two excess dispersions
  differ. To impose $`\phi_{x1} = \phi_{x2}`$ instead, route both
  through one non-linear parameter — `nlf(shapexone ~ shapexx)`,
  `nlf(shapextwo ~ shapexx)`, `shapexx ~ 1`, with `nl = TRUE` — and set
  its prior with `nlpar = "shapexx"`. That is the model this package
  fitted under a single `shapex` dpar before 0.10.0, term for term.
- [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  — for data where the first count is missing on some rows. It is the
  same family: same name, same six dpars, same likelihood and the same
  post-processing. What it adds is a second `vint()` integer, a 0/1 flag
  saying whether $`y_1`$ was recorded. Rows where it was not are scored
  by the integrated-out marginal of this same joint rather than dropped,
  so they still inform $`\mu`$, $`\phi_s`$, $`\lambda_2`$, $`\phi_{x2}`$
  and any group-level effects; afterwards the fit imputes the missing
  count conditional on the observed one. Note that $`\lambda_1`$ and
  $`\phi_{x1}`$ are then identified by the matched rows alone. Despite
  its former name, this is unrelated to brms’s `cens()` addition term,
  which means a value known to lie in a set.
- The $`(M, f, \delta)`$**reparameterisation** — overall level $`M`$,
  congruence $`f`$, and source bias $`\delta`$, with the dispersions on
  an SD scale $`\kappa = 1/\sqrt{\phi}`$ where $`\kappa = 0`$ is the
  Poisson limit. It is fitted through
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  rather than a separate family, and
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md)
  converts between the two. Because $`\delta = 0`$ and $`\kappa = 0`$
  are finite, interpretable nulls, these coordinates admit shrinkage
  priors towards the simpler model unless the data support otherwise
  ([Simpson et al. 2017](#ref-simpsonPenalisingModelComponent2017)). The
  rate parameterisation has no such nulls. Worked end to end in the
  article [*The anatomy of a paired
  count*](https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.html).
- [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  — the same partial observation, equidispersed:
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)’s
  Poisson special case. Its integrated-out marginal is closed form,
  $`y_2 \sim \mathrm{Poisson}(\mu + \lambda_2)`$, since a sum of
  independent Poissons is Poisson. Use it when the margins are not
  overdispersed, rather than fitting
  [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  with its dispersions pressed against the Poisson boundary.
- **Modelling the difference instead.** Where only the disagreement
  $`d = y_1 - y_2`$ is of interest, or where truncation via
  [`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html)
  is required, the Skellam and discrete Laplace/normal families in the
  companion package
  [`skellambrms`](https://github.com/anhsmith/skellambrms) model $`d`$
  directly. The difference of two independent Poisson counts is
  Skellam-distributed ([Skellam 1946](#ref-skellam1946)); for the
  Bayesian treatment of count differences generally, see Karlis and
  Ntzoufras ([2006](#ref-karlis2006)). Note what that reduction
  discards: the shared component cancels out of $`d`$, so the level and
  the congruence are no longer estimable.

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

Kirkpatrick, Robert M. 2022. *RMKdiscrete: Sundry Discrete Probability
Distributions*. <https://CRAN.R-project.org/package=RMKdiscrete>.

Kirkpatrick, Robert M., and Michael C. Neale. 2016. “Applying
Multivariate Discrete Distributions to Genetically Informative Count
Data.” *Behavior Genetics* 46 (2): 252–68.
<https://doi.org/10.1007/s10519-015-9757-z>.

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
