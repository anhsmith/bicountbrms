# bicountbrms <img src="man/figures/logo.png" align="right" height="139" alt="bicountbrms logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/anhsmith/bicountbrms/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Custom [brms](https://paulbuerkner.com/brms/) families (Bürkner 2017) for
modelling a **pair of counts** from two sources that are meant to measure the
same thing — two observers, two instruments, two reporting channels — where the
question is how much, and how systematically, they disagree.

The pair $(y_1, y_2)$ is modelled itself, not reduced to its difference, using
**trivariate reduction** (Holgate 1964): the two counts share an unobserved
latent component, which is what makes them correlated, plus a private component
each, whose difference is the disagreement. One model then answers:

- **how much was there?** — the pair's overall level, not just its difference;
- **how much did the two sources agree?** — their congruence, how much of
  the total each saw in common;
- **was there bias toward one source?**, and by how much;
- **how overdispersed are the counts?** — excess dispersion relative to a
  Poisson;
- **what would the other source have recorded?** — often one source is complete
  and the other is sometimes not observed; the model simulates $y_1$ given an
  observed $y_2$, including on rows where $y_1$ was never observed at all.

Two families — equidispersed and overdispersed — each with two constructors,
according to whether the first count was recorded on every row:

| | both counts observed | first count sometimes unobserved |
|---|---|---|
| **equidispersed** (Poisson latents) | `bipois()` | `bipois_partialobs()` |
| **overdispersed** (negative-binomial latents) | `binegbin()` | `binegbin_partialobs()` |

Each row of that table is one family: both constructors return the same
likelihood and the same post-processing methods, and differ only in whether the
data carries an observation flag. A fully paired design never has to supply one.

The three rates the likelihood takes can be reparameterised as overall level
$M$, congruence $f$, and source bias $\delta$, which separates the first three
questions above. Predictors and random effects can be placed on each
independently, and $\delta = 0$ and $\kappa = 0$ are finite nulls. See
[The $(M, f, \delta)$ coordinates](#the-m-f-delta-coordinates) below, and the
article [*The anatomy of a paired count*][anatomy] for a worked fit.

**Modelling the difference instead.** Where only the disagreement
$d = y_1 - y_2 \in \mathbb{Z}$ is of interest, or where truncation via
`resp_trunc()` is required, the companion package
[`skellambrms`](https://github.com/anhsmith/skellambrms) supplies Skellam,
discrete Laplace and discrete normal families that model $d$ directly. Both
packages were formerly distributed together as `pairedcountbrms`; see
[Installation](#installation).

## Overview

Two counts of the same quantity are correlated, and no bivariate count
distribution parameterises that correlation directly, as the bivariate normal
does. Dependence must therefore be constructed. Under trivariate reduction it
follows from the counting process itself: the two sources record overlapping
sets of events, and the dependence is the count they share. This admits only
positive correlation. Copulas are an alternative, imposing a dependence
structure on transformed margins rather than generating one from the count
parameters (Genest and Nešlehová 2007).

Modelling each count separately assumes that dependence is zero. Standard
`brms` count families (Poisson, negative binomial, …) model a single
non-negative response, and cannot take a jointly-modelled pair.

This package supplies two families with likelihood $P(y_1, y_2)$.

The `_partialobs()` constructors admit rows on which $y_1$ was not observed. On
those rows the likelihood is the same joint with $y_1$ integrated out over its
whole support, so they still inform the shared component and the second
source's rate rather than being dropped. After fitting, `posterior_predict()`
returns a posterior distribution for $y_1$ on every row, conditional on that
row's observed $y_2$ — including the rows where $y_1$ was never recorded.

This is **not** censoring in the sense of brms's `cens()` addition term, which
means a value known to lie in a set ($\le c$, $\ge c$, or within an interval).
Here $y_1$ is not observed at all. These constructors were named `_cens` up to
0.9.1 and renamed at 0.10.0 for that reason; do not combine them with `cens()`.

For the bivariate Poisson the difference $y_1 - y_2$ is exactly
Skellam-distributed (Skellam 1946), so `bipois()` induces the Skellam as its
difference model.

## The construction

Both families are built by the same **trivariate reduction** (Holgate 1964;
Karlis and Ntzoufras 2003): three independent latent counts

$$
N_{\text{shared}} \sim \mathcal{D}(\mu), \quad
N_1 \sim \mathcal{D}(\lambda_1), \quad
N_2 \sim \mathcal{D}(\lambda_2),
$$

combined as

$$
y_1 = N_{\text{shared}} + N_1, \qquad
y_2 = N_{\text{shared}} + N_2.
$$

The shared count $N_{\text{shared}}$ appears in both, inducing positive
correlation; the two private counts $N_1, N_2$ drive the difference. Source 1 is
the modelled response; source 2 is supplied via `vint()`. $\mu$ is the shared
*rate*, **not** the mean of either response. The likelihood marginalises the
unobserved $N_{\text{shared}}$ analytically:

$$
P(y_1=x,\, y_2=y)
  = \sum_{k=0}^{\min(x,y)}
    f_{\text{s}}(k)\, f_1(x-k)\, f_2(y-k),
$$

where $f_{\text{s}}, f_1, f_2$ are the pmfs of the three latent counts. Two
consequences follow directly:

- **The margins.** $\mathrm{E}[y_1] = \mu + \lambda_1$ and likewise for $y_2$;
  $\mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_{\text{shared}})$, so the correlation
  is
  $\mathrm{Var}(N_{\text{shared}}) / \sqrt{\mathrm{Var}(y_1)\mathrm{Var}(y_2)}$.
- **The difference.** $d = y_1 - y_2 = N_1 - N_2$ — the shared count *cancels*.
  So the difference depends only on the two private components. For the Poisson
  case that difference is precisely $\mathrm{Skellam}(\lambda_1, \lambda_2)$.

### The $(M, f, \delta)$ coordinates

The three rates enter the likelihood directly, but none of them is separately
interpretable: raising the overall level of counting moves all three together.
Three coordinates separate the quantities of interest:

| | | |
|---|---|---|
| $M$ | overall level | $\mu + (\lambda_1 + \lambda_2)/2$ |
| $f$ | congruence — the share of $M$ both sources saw | $\mu / M$ |
| $\delta$ | source bias, on a log-ratio scale | $\tfrac12\log(\lambda_1/\lambda_2)$ |

with dispersions on an SD scale, $\kappa = 1/\sqrt{\phi}$, so that $\kappa = 0$
is the Poisson limit. [`binegbin_mfd_to_dpars()`][mfd] and
[`binegbin_dpars_to_mfd()`][mfd] convert both ways; the map is a bijection. The
rate half applies to both families and both constructors; only the dispersion
half is specific to the negative-binomial one.

These coordinates separate effects that the rates confound. A group-level effect
on $\mu$ changes both the level and the congruence; an effect on $M$ or on $f$
changes one at a time. They also supply finite nulls — $\delta = 0$ is no bias,
$\kappa = 0$ is Poisson — and so admit **penalised-complexity priors** (Simpson
et al. 2017), which shrink to the simpler model unless the data support
otherwise. The rate parameterisation has no such nulls: its Poisson limit is
$\phi \to \infty$.

No separate family is required; the coordinates are supplied through `nlf()`.
The article [*The anatomy of a paired count*][anatomy] simulates from known
$(M, f, \delta)$, fits, recovers the coordinates, and converts back, alongside
an interactive widget linking the two parameterisations.

[mfd]: https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.html
[anatomy]: https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.html

### Poisson or negative-binomial latents: overdispersion

| Latent law | Family | Dispersion dpars | Var of each latent | Use when |
|---|---|---|---|---|
| Poisson | `bipois()` | none | $\mathrm{Var}=\text{mean}$ | margins are **not** overdispersed |
| negative-binomial | `binegbin()` | `shapes` (shared), `shapexone`/`shapextwo` (private) | $m + m^2/\phi$ | margins **are** overdispersed (the usual case) |

With Poisson latents, each component has $\mathrm{Var}=\text{mean}$, so
`bipois()` cannot represent overdispersed margins and underfits the marginal
(and difference) variance of real count data. `binegbin()` replaces each latent
with $N \sim \text{NB2}(m, \phi)$ — Stan's `neg_binomial_2`, R's
`dnbinom(size = φ, mu = m)` — with mean $m$ and variance $m + m^2/\phi$. It
carries the extra spread in **scalar** dispersion dpars: `shapes`
$=\phi_{\text{s}}$ for the shared count, and `shapexone`/`shapextwo`
$=\phi_{\text{x}1},\phi_{\text{x}2}$, one per private count. The moments become

$$
\mathrm{Var}(y_1) = \Big(\mu+\tfrac{\mu^2}{\phi_{\text{s}}}\Big) +
  \Big(\lambda_1+\tfrac{\lambda_1^2}{\phi_{\text{x}1}}\Big),
\qquad
\mathrm{Cov}(y_1,y_2) = \mu+\tfrac{\mu^2}{\phi_{\text{s}}},
$$

and, with $\lambda_1=\lambda_2=\lambda$ and
$\phi_{\text{x}1}=\phi_{\text{x}2}=\phi_{\text{x}}$,
$\mathrm{Var}(d) = 2\big(\lambda + \lambda^2/\phi_{\text{x}}\big)$. As every
$\phi\to\infty$ the negative-binomials collapse to Poissons and `binegbin`
$\to$ `bipois`, under either constructor. That limit is not only a description:
it is used as a test (see [Testing](#testing)).

**Two private dispersions, and how to tie them.** The two sources are different
instruments, and there is no reason their source-only excess must be equally
overdispersed, so `binegbin()` estimates one for each. Before 0.10.0 it carried
a single `shapex` and only the partially observed family had the pair, which
was the wrong way round: every row of a fully paired design informs both, so
that is the case in which they are easiest to identify. To impose
$\phi_{\text{x}1}=\phi_{\text{x}2}$, route both through one non-linear
parameter and set the prior with `nlpar = "shapexx"`:

``` r
bf(y1 | vint(y2) ~ 1,
   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
   shapexx ~ 1, ..., nl = TRUE)
```

Under **partial** observation the two are not equally well identified:
$\phi_{\text{x}2}$ governs the always-observed margin and appears on both
branches, while $\phi_{\text{x}1}$ appears only on the matched branch and is
identified by the matched rows alone.

Where the counts really are equidispersed, prefer the Poisson family to a
negative-binomial one fitted with its dispersions pressed against
$\phi \to \infty$, which is a boundary at which sampling degrades. Compare the
two with `loo()`.

**Why scalar dispersion, not a random effect.** An observation-level random
effect on the private components fails synthetic recovery: with one observed
pair but three latent deviates per unit, the excess deviates absorb residual
variation and fresh draws do not regenerate the observed spread. A conditional
posterior-predictive check does not reveal this; a marginal one does. Scalar
`shapes`/`shapexone`/`shapextwo` are instead identified from the aggregate mean–variance
relationship across units. The `binegbin.R` file header gives the detail.

### `bipois_partialobs()` and `binegbin_partialobs()`: partially observed pairs

These constructors admit rows on which $y_1$ was **not observed at all**. Each
row carries a second `vint()` integer, a `y1_obs` $\in\{0,1\}$ flag — `1` where
both counts were recorded, `0` where the first was not. Put any non-negative
integer in `y1` on those rows (`0` is the conventional placeholder); the
likelihood does not read it. Do not use `NA`, which brms drops before fitting,
taking the row's observed $y_2$ with it.

- **`y1_obs == 1` (matched):** the full joint lpmf on $(y_1, y_2)$ — the same
  function `bipois()`/`binegbin()` calls, not a second copy of it.
- **`y1_obs == 0` ($y_2$-only):** the $y_1$-**integrated marginal** of that
  *same* joint,

$$
P(y_2=y) = \sum_{k=0}^{y} f_{\text{s}}(k)\, f_2(y-k),
$$

  i.e. the joint with $N_1$ summed out (the inner sum over $y_1$ telescopes to
  $1$). This is a convolution of the shared and source-2 private components —
  **not** a separate univariate family on $y_2$, which would be a different,
  incoherent model.

One `brm()` call thus pools matched and unmatched rows under one likelihood. An
unmatched row is not dropped and is not given a different model: it still
informs the shared component, the second source's rate and dispersion, and any
group-level effects. What it cannot inform is $\lambda_1$, `shapexone` and the
between-source bias, which appear only on the matched branch and are identified
by the matched rows **alone**. A design with 20 matched rows in 500 will learn
them weakly and lean on their priors — worth knowing before collecting data,
not after seeing a posterior. See [Limitations](#limitations).

After fitting, `posterior_predict()` and `posterior_epred()` return a draw and
$\mathrm{E}[y_1 \mid y_2]$ for **every** row, matched and unmatched alike:
`y1_obs` selects a likelihood branch, not a prediction. Imputing the unobserved
first count conditional on the observed second one is usually the reason for
wanting this.

**Which shape does a stored fit have?** `family$name` is the same either way, so
the flag's presence is what distinguishes them:

```r
"vint2" %in% names(standata(fit))   # TRUE for a partially observed fit
fit$family$vars                     # c("vint1[n]", "vint2[n]") or c("vint1[n]", "1")
```

**The Poisson case's marginal is closed form.** A sum of independent Poissons is
Poisson, so for `bipois_partialobs` the convolution above collapses exactly to

$$
y_2 \sim \mathrm{Poisson}(\mu + \lambda_2),
$$

with no sum to evaluate. `binegbin_partialobs` has no such shortcut, because
$\mathrm{NB2} + \mathrm{NB2}$ is not $\mathrm{NB2}$ unless the two components
share a success probability, so it evaluates the convolution term by term. This
also makes `bipois_partialobs` the analytic reference against which
`binegbin_partialobs`'s unmatched branch is checked.

It sharpens what the unmatched rows can say, too. For `bipois_partialobs` they see
$\mu$ and $\lambda_2$ only through their **sum**: they constrain the total rate
of the observed margin, not how it splits between shared and private. Separating
the two — and so estimating the congruence $f$ — rests on the matched rows
however many unmatched rows there are.

**The two excess dispersions under partial observation.** `binegbin()` and
`binegbin_partialobs()` both give each private component its own dispersion —
`shapexone` $=\phi_{\text{x}1}$ for $N_1$ and `shapextwo` $=\phi_{\text{x}2}$
for $N_2$ — six dpars in all. The moments are the same either way:

$$
\mathrm{Var}(y_1) = \Big(\mu+\tfrac{\mu^2}{\phi_{\text{s}}}\Big) +
  \Big(\lambda_1+\tfrac{\lambda_1^2}{\phi_{\text{x}1}}\Big),
\qquad
\mathrm{Var}(y_2) = \Big(\mu+\tfrac{\mu^2}{\phi_{\text{s}}}\Big) +
  \Big(\lambda_2+\tfrac{\lambda_2^2}{\phi_{\text{x}2}}\Big),
$$

while $\mathrm{Cov}(y_1,y_2) = \mu+\mu^2/\phi_{\text{s}}$ is unchanged, since
only the shared component contributes to it. Freeing the two dispersions
therefore alters the margins and the difference but not the covariance, and the
correlation changes only through its denominator. For the difference,
$\mathrm{Var}(d) = \big(\lambda_1+\lambda_1^2/\phi_{\text{x}1}\big) +
\big(\lambda_2+\lambda_2^2/\phi_{\text{x}2}\big)$, which at
$\lambda_1=\lambda_2=\lambda$ reduces to
$2\lambda + \lambda^2\big(1/\phi_{\text{x}1} + 1/\phi_{\text{x}2}\big)$: the two
dispersions enter through their reciprocals, so the more overdispersed source
contributes the greater share of the disagreement. Setting
$\phi_{\text{x}1}=\phi_{\text{x}2}$ recovers the single-dispersion expressions
of earlier versions.

What differs between the two constructors is not the parameters but how well
they are identified. `shapextwo` governs the always-observed margin and so
appears on both branches, and the unmatched rows sharpen it. `shapexone`
appears only on the matched branch, exactly as $\lambda_1$ does. Under full
pairing that asymmetry vanishes, which is why `binegbin()` — not
`binegbin_partialobs()` — is the constructor for which two free dispersions are
easiest to estimate.

### Usage

```r
library(brms)
library(bicountbrms)

# bipois(): joint bivariate Poisson (non-overdispersed margins)
fit_bp <- brm(
  bf(y1 | vint(y2) ~ 1,
     mu ~ 1 + (1 | vessel), lambdaone ~ 1, lambdatwo ~ 1),
  data = dat, family = bipois(), stanvars = bipois_stanvars(), chains = 4
)

# binegbin(): joint bivariate negative-binomial (overdispersed margins)
fit_nb <- brm(
  bf(y1 | vint(y2) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx), lamx ~ 1,
     shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE),
  data = dat, family = binegbin(), stanvars = binegbin_stanvars(), chains = 4
)

# ...or its symmetric special case, one dispersion for both private components
fit_nb_sym <- brm(
  bf(y1 | vint(y2) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx), lamx ~ 1,
     nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
     shapes ~ 1, shapexx ~ 1, nl = TRUE),
  data = dat, family = binegbin(), stanvars = binegbin_stanvars(), chains = 4
)

# bipois_partialobs(): equidispersed, and y1 is unobserved where y1_obs == 0
fit_pj <- brm(
  bf(y1 | vint(y2, y1_obs) ~ 1,
     mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
     nlf(lambdaone ~ lamx + methd),
     nlf(lambdatwo ~ lamx - methd),
     lamx ~ 1, methd ~ 1, nl = TRUE),
  data = dat, family = bipois_partialobs(), stanvars = bipois_partialobs_stanvars(),
  chains = 4
)

# binegbin_partialobs(): the same partial observation, overdispersed margins
fit_cj <- brm(
  bf(y1 | vint(y2, y1_obs) ~ 1,
     mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
     nlf(lambdaone ~ lamx + methd),
     nlf(lambdatwo ~ lamx - methd),
     lamx ~ 1, methd ~ 1, shapes ~ 1,
     shapexone ~ 1, shapextwo ~ 1, nl = TRUE),
  data = dat, family = binegbin_partialobs(), stanvars = binegbin_partialobs_stanvars(),
  chains = 4
)

# ...or its symmetric special case, one dispersion for both private components
fit_cj_sym <- brm(
  bf(y1 | vint(y2, y1_obs) ~ 1,
     mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
     nlf(lambdaone ~ lamx + methd),
     nlf(lambdatwo ~ lamx - methd),
     nlf(shapexone ~ shapexx),
     nlf(shapextwo ~ shapexx),
     lamx ~ 1, methd ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE),
  data = dat, family = binegbin_partialobs(), stanvars = binegbin_partialobs_stanvars(),
  chains = 4
)
```

The `nlf(lambdaone ~ lamx)` / `nlf(lambdatwo ~ lamx)` idiom ties the two private
rates to one value — a "no systematic bias" assumption,
$\mathrm{E}[y_1]=\mathrm{E}[y_2]$. Splitting them as `lamx + methd` /
`lamx - methd` introduces a directional bias parameter `methd` (half the log
ratio of the two private rates); giving the two rates separate predictors
(`nlf(lambdaone ~ l1)`, `nlf(lambdatwo ~ l2)`, `l1 ~ 1`, `l2 ~ 1`) is the fully
unconstrained version.

The second count (and, for the `_partialobs()` constructors, the flag) travels
via `vint()` because `custom_family()` declares a single response column —
`vint(y2, y1_obs)` binds `vint1 = y2`, `vint2 = y1_obs` in listed order. On
unmatched rows the response column is not read by the likelihood, so any
integer placeholder there is inert; supply one rather than `NA`, which `brms`
drops before the family is reached, taking the observed $y_2$ with it.

**A fully paired model never mentions the flag.** The plain constructors declare
`vars = c("vint1[n]", "1")`, and `brms` pastes those entries into the generated
Stan call verbatim, so `binegbin()` reaches the same `binegbin_lpmf` with
`y1_obs` fixed at `1`:

```
binegbin()             target += binegbin_lpmf(Y[n] | ..., vint1[n], 1);
binegbin_partialobs()  target += binegbin_lpmf(Y[n] | ..., vint1[n], vint2[n]);
```

That is why there is one Stan function per family rather than two overloaded
ones, and so why the package declares no floor on the Stan version — function
overloading needs Stan 2.29.

### Set priors on the dispersions

`brms` assigns a custom family's **`mu` dpar** a default `student_t` prior but
leaves the remaining dpars **flat and improper**; `get_prior()` returns an empty
`prior` column for `lamx`, `shapes`, `shapexone` and `shapextwo`. These are the
weakly-identified parameters, and unregularised the chains explore the flat
tail. In this package's recovery test, omitting the priors produced a divergent
transition and $\hat R = 1.0101$, although all parameters still recovered.
Weakly informative priors are sufficient:

```r
prior = c(
  prior(normal(2, 1), class = "Intercept"),                      # log shared rate
  prior(normal(2, 1), class = "b",         nlpar = "lamx"),      # log private rate
  prior(normal(2, 1), class = "Intercept", dpar  = "shapes"),
  prior(normal(2, 1), class = "Intercept", dpar  = "shapexone"),
  prior(normal(2, 1), class = "Intercept", dpar  = "shapextwo")
)
```

Under the symmetric `nlf(shapexone ~ shapexx)` idiom the last two lines become
`prior(normal(2, 1), class = "b", nlpar = "shapexx")` instead — routing a dpar
through a non-linear parameter moves its prior from `dpar =` to `nlpar =`.
`bipois()` and `bipois_partialobs()` have no dispersion dpars, so only the rate
lines apply to them.

These suit rates and dispersions of roughly 1 to 50. To adapt them, shift the
**mean** to the scale of the counts rather than increasing the SD. Every dpar of
every family here is log-linked, so a normal prior is lognormal on the natural
scale, and increasing its SD moves mass to implausible values rather than making
the prior neutral (Smith et al. 2020, supplement 3 — see
[References](#references)).

The `class`/`dpar`/`nlpar` slots differ between a rate supplied through `nlf()`
(`class = "b"`, `nlpar`) and the dispersions (`class = "Intercept"`, `dpar`).
`get_prior()` returns the required combination.

Under the $(M, f, \delta)$ parameterisation the dispersions admit stronger
justification. On the SD scale $\kappa = 1/\sqrt{\phi}$ the Poisson limit is a
finite $\kappa = 0$, so `exponential()` is a penalised-complexity prior that
shrinks to Poisson unless the data support overdispersion, and
`double_exponential()` on $\delta$ shrinks to no between-source bias.
[*The anatomy of a paired count*][anatomy] demonstrates both.

### Prediction and expectation

`posterior_predict()` and `posterior_epred()` answer the same question — what
the first source would have recorded, **given the second count that was actually
observed** — one by simulation, the other in expectation. $y_2$ is fixed data
and is never itself re-simulated.

$$
\mathrm{E}[y_1 \mid y_2] = \mathrm{E}[N_{\text{shared}} \mid y_2] + \lambda_1
$$

- **`bipois` / `bipois_partialobs`:** the conditional split is closed form.
  Conditioning a sum of independent Poissons on its total gives a Binomial, so
  $N_{\text{shared}}\mid y_2 \sim \text{Binomial}\big(y_2,\,
  \mu/(\mu+\lambda_2)\big)$, then a fresh $N_1$; the expectation is
  $y_2\,\mu/(\mu+\lambda_2) + \lambda_1$.
- **`binegbin` / `binegbin_partialobs`:** a negative-binomial sum has no Binomial
  conditional, so the discrete law $P(N_{\text{shared}}=k\mid y_2)
  \propto f_{\text{s}}(k)\,f_2(y_2-k)$ over $k = 0,\ldots,y_2$ is sampled
  directly for prediction, and summed for the expectation. Both are exact.
- Under `_partialobs()`, `y1_obs` is **ignored** by both: every row (matched
  and unmatched alike) gets a conditional $y_1$, so the unobserved margin can be
  imputed fleet-wide. That is the point of fitting them.

Every family's `posterior_epred_<family>()` agrees with the mean of its own
`posterior_predict_<family>()` draws to within Monte Carlo error, which is the
standard the test suite holds them to.

Both are called the ordinary way — `posterior_epred(fit)`, and hence `fitted()`
and `conditional_effects()`, dispatch to the family's own method, which a
package test verifies against a direct
`posterior_epred_<family>(prepare_predictions(fit))` call. (`brms` checks
truncation before family type and its truncated branch has no custom-family
fallback, so `posterior_epred()` errors on a truncated custom family. That
cannot arise here, since these families do not support `resp_trunc()`; it does
affect `skellambrms`, which documents the workaround.)

## Parameterisation and naming notes

**The forced `"mu"` dpar.** `brms::custom_family()` unconditionally requires one
dpar to be named literally `"mu"`. Here that slot holds the shared *rate*, not a
response mean: $\mathrm{E}[y_1] = \mu + \lambda_1$, so `mu` is neither the mean
of either count nor the mean of their difference. If you read `make_stancode()`
output or call `get_dpar(prep, "mu")`, that is what you are looking at.

**`lambdaone`, not `lambda1`.** `custom_family()` rejects any dpar name ending
in a digit (`"'dpars' should not end with a number."`), as well as dots and
underscores. The two private rates are therefore spelled `lambdaone` and
`lambdatwo` in code, while the documentation writes them $\lambda_1$ and
$\lambda_2$. Nothing else distinguishes the two spellings — they are the same
parameter.

**Under `nl = TRUE`, give `mu` its own formula.** The examples above use
`nl = TRUE` so that `nlf()` can tie the two private rates together. In a
non-linear `brms` formula the main formula's right-hand side is a *non-linear
expression* for `mu`, **not** a request for an intercept. So this:

```r
bf(y1 | vint(y2) ~ 1,                      # <- no `mu ~ 1` term
   nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx),
   lamx ~ 1, shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE)
```

generates `mu[n] = exp(1);` — the shared rate is **pinned at $e \approx 2.72$ and
never estimated**. It does not error. `brms` simply fits a different model, and
because the pair's level then has to be absorbed by the private rates, sampling
becomes very slow. Add `mu ~ 1` (or `mu ~ 1 + (1 | group)`, as in the examples
above) and the generated code becomes `mu += Intercept; mu = exp(mu);` as
intended. Checking `make_stancode()` for a `b_Intercept` is the quickest way to
confirm.

### Notation

| Code | Math | Meaning |
|---|---|---|
| `y1` (the response) | $y_1$ | count from source 1 |
| `y2` (via `vint()`) | $y_2$ | count from source 2 |
| — | $N_{\text{shared}}$ | latent count both sources saw |
| — | $N_1$, $N_2$ | latent counts only one source saw |
| `mu` | $\mu$ | rate of the shared component (**not** a response mean) |
| `lambdaone`, `lambdatwo` | $\lambda_1$, $\lambda_2$ | rates of the two private components |
| `shapes` | $\phi_{\text{s}}$ | NB2 dispersion of the shared component |

| `shapexone`, `shapextwo` | $\phi_{\text{x}1}$, $\phi_{\text{x}2}$ | NB2 dispersions of the two private components (`binegbin`) |
| `y1_obs` (via `vint()`) | — | the `_partialobs()` flag: was $y_1$ observed? |

Source 1 is whichever count you put on the left-hand side of the formula; the
labelling carries no further meaning.

## Installation

```r
# install.packages("pak")
pak::pak("anhsmith/bicountbrms")
```

Stan and a C++ toolchain are required. On Windows, install
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).
Works with either rstan or cmdstanr as the brms backend.

**Package history.** These families were released as `skellambrms` (0.1.0–0.5.0)
and then `pairedcountbrms` (0.6.0–0.8.0), which carried both these joint
families and a set of families for modelling the difference $d = y_1 - y_2$
directly. At 0.9.0 the two were separated: the joint families are this package,
and the difference families returned to
[`skellambrms`](https://github.com/anhsmith/skellambrms). The split itself
renamed nothing, so `library(pairedcountbrms)` becomes `library(bicountbrms)` or
`library(skellambrms)` and that is all it requires.

`github.com/anhsmith/pairedcountbrms` redirects here.

**Migrating from 0.9.1 or earlier.** 0.10.0 is a breaking release; `NEWS.md`
carries the full account.

The family surface changed. `binegbin_joint()` became `binegbin_cens()` at
0.9.0; at 0.10.0 `bipois_cens()` and `binegbin_cens()` became
`bipois_partialobs()` and `binegbin_partialobs()`, and all seven old names were
**removed**, with no deprecation layer. In source code, substitute the new name.

| up to 0.9.1 | from 0.10.0 |
|---|---|
| `bipois_cens()` | `bipois_partialobs()` |
| `binegbin_cens()`, `binegbin_joint()` | `binegbin_partialobs()` |
| `binegbin()` with `shapex` | `binegbin()` with `shapexone`/`shapextwo` |

`binegbin()` gained the second excess dispersion, so `shapex ~ 1` becomes
`shapexone ~ 1, shapextwo ~ 1`, or the tied
`nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1` for the
old single-dispersion model. A prior set with `dpar = "shapex"` moves to
`nlpar = "shapexx"` under the tied form.

**Fits stored under the old names need the old package.** A fit records its
`family$name` permanently and brms builds `log_lik_<name>` from that string at
call time, so a fit made with `binegbin_cens()` or `binegbin_joint()` needs a
version that still defines those methods. Install 0.9.1 and keep it attached
for that fit; nothing you pass to `loo()` or `posterior_predict()` can route
around it, and 0.10.0 deliberately ships no forwarders.

Stored **`binegbin`** fits are the exception and need nothing: they carry the
name this release keeps, and their single `shapex` resolves to both `shapexone`
and `shapextwo` on the current methods, as pre-0.7.0 `lambdaem`/`lambdalb`
resolve to `lambdaone`/`lambdatwo`. A package test pins that exact combination.

One numerical result also changed, at 0.9.0: see `posterior_epred_binegbin()`
under [Limitations](#limitations), and `NEWS.md`.

Documentation, including a worked getting-started vignette that simulates,
fits, and recovers `binegbin` parameters end to end, is at
<https://anhsmith.github.io/bicountbrms/>. Locally:

```r
vignette("bicountbrms")
```

## Limitations

**The correlation cannot be negative.** It is generated by a shared component
contributing non-negatively to both counts, so a design in which one source
systematically records *more* when the other records less — competing rather
than concurrent observation — is outside what these families can represent.

**No absolute standard is identified.** $M$ is a midpoint of what the two
sources *report*. Events both sources miss are invisible to the model, and no
amount of data recovers them; a claim about the true level requires an external
reference the model does not supply.

**Congruence and co-detection are different quantities.** $f$ is
mass-weighted — the share of expected count both sources recorded together —
whereas the proportion of units on which both recorded anything counts units.
Adding one-sided mass on units where both already record moves the two in
opposite directions, so they should not be substituted for one another.

**Truncation is not supported.** `resp_trunc()` does not apply to these
families; no `_lccdf_stanvars()` is provided for them. The difference families
in [`skellambrms`](https://github.com/anhsmith/skellambrms) do support it.

**Identifiability under sparse pairing.** $\lambda_1$ and
the between-source bias are informed only by the matched (`y1_obs == 1`) rows —
as is `shapexone` for `binegbin_partialobs` — so a fit with few matched rows will
learn them weakly even if the total sample is large. The unmatched rows add power
for the shared and level parameters and for `shapextwo`, not for the bias. For
`bipois_partialobs` the unmatched rows see $\mu$ and $\lambda_2$ only through their
sum, so the congruence $f$ is likewise a matched-row quantity. Fits that lean on
the bias, on the congruence, or on a difference between the two excess
dispersions should be judged against the matched subset, not the full $n$.

**`posterior_epred_binegbin()` changed at 0.9.0.** Earlier versions substituted
the *marginal* shared fraction $\mu/(\mu+\lambda_2)$ for the *conditional* one —
the `bipois` answer, exact only in the Poisson limit. Over a grid of plausible
rates and dispersions, with $y_2$ within one SD of its mean, that substitution
was in error by more than 5% in about half the settings and more than 20% in a
third, in either direction depending on which component carried the greater
dispersion. It is now the exact conditional expectation. Numbers computed with
an earlier version's `posterior_epred_binegbin()` should be recomputed;
`log_lik()`, `loo()` and `posterior_predict()` are unaffected.

## Testing

For both families the suite (`tests/testthat/`) validates:

- the marginalised joint log-PMF against an **independent R brute-force
  reference** across a rate/shape grid and at edge cases, via
  `rstan::expose_stan_functions()` (agreement to ~$10^{-14}$);
- **normalisation** to $1$, over the pair for the matched branch and over $y_2$
  for the unmatched one;
- the analytic **moment identities** — mean, marginal variance, difference
  variance, covariance;
- the **Poisson-limit reduction** `binegbin` $\to$ `bipois`, on both
  branches, checked as a limit (the
  error must shrink as $\phi$ grows) rather than at a single tolerance;
- for `bipois`, that the induced difference $d = y_1 - y_2$ matches
  `skellam::dskellam()`;
- end-to-end **parameter recovery** with divergence/$\hat R$ checks.

For the partially observed branch it additionally pins the **marginal identity**
($\sum_{y_1}$ of the matched branch equals the $y_2$-only branch), the
agreement of the Stan lpmf with the R reference at both values of the flag,
the **conditional-prediction identity** (`posterior_predict` draws
match joint/marginal) on matched **and unmatched** rows and with the two excess
dispersions distinct, and — for `bipois_partialobs` — that the closed-form
unmatched branch equals the brute-force convolution, which is what licenses
using the closed form at all.

`test-epred.R` states the `posterior_epred` convention once for both
families and holds each to the same standard: the expectation must equal the
mean of that family's own `posterior_predict` draws.

`test-partialobs-predict.R` covers the prediction path on **unmatched**
rows with `shapexone` and `shapextwo` distinct — the combination under which
exchanging the two inside `posterior_predict_binegbin()` passed every other
test in the suite. Its closing comment records what is still not covered, and
why one candidate check is deliberately absent rather than overlooked.

### Running the tests

Everything that needs a Stan toolchain is behind `skip_on_cran()`, and every
file-level compilation behind `stan_tests_enabled()`, so a plain `R CMD check`
runs the analytic and R-side tests but **skips every model fit and every
compile**. To include them:

```bash
NOT_CRAN=true R CMD check --no-manual bicountbrms_0.10.0.tar.gz
```

The recovery tests distinguish two questions.

- **Smoke gates** run with the other fitting tests. Each fits once and checks
  convergence (zero divergences, $\hat R < 1.02$) and that the true value lies
  within a wide (99%) posterior interval, detecting gross mis-specification.
  This is not a calibration statement: for a correct model with a calibrated
  posterior, asserting that a truth falls within a 90% interval from a single
  fit fails 10% of the time by construction, and across ~18 such assertions
  spurious failures were the norm.
- **Calibration** is assessed separately, as the proportion of repeated
  simulate-and-refit replicates whose nominal interval contains the truth. This
  costs minutes, and is opt-in:

```bash
BICOUNTBRMS_COVERAGE=true NOT_CRAN=true Rscript -e 'testthat::test_local()'
```

The pass threshold is derived from the $\mathrm{Binomial}(R, 0.9)$ null rather
than chosen by eye, so the false-failure rate is a stated quantity: 0.16% at
$R = 10$. Run it before a release, or after changing a likelihood, link, or
parameterisation. See `tests/testthat/helper-coverage.R`.

**What runs automatically.** `.github/workflows/R-CMD-check.yaml` defines two
jobs on different cadences, since the checks differ by an order of magnitude in
cost:

| job | when | what |
|---|---|---|
| `check` | every push and pull request | builds, docs, examples, and the tests that need no Stan backend |
| `check-stan` | weekly, and on demand | sets `NOT_CRAN=true` and installs both backends, so every model fit actually runs |

The fitting tests are both the most likely to detect a regression and too slow
for every push. Leaving them behind `skip_on_cran()` with nothing scheduled to
lift it means they do not run at all: two assertions in this package remained
broken undetected for that reason.

Trigger `check-stan` by hand from the repository's Actions tab whenever you want
it. GitHub disables scheduled workflows after 60 days without a commit, so if
this package goes quiet the weekly run stops; the manual trigger is the fallback.

## Function reference

Each family exports two constructors, their `_stanvars()`, and one set of
`log_lik_`, `posterior_predict_` and `posterior_epred_` interface functions.

| Function | Purpose |
|---|---|
| `bipois()` / `bipois_stanvars()` | Joint bivariate Poisson, both counts observed |
| `bipois_partialobs()` / `bipois_partialobs_stanvars()` | The same, admitting rows where $y_1$ was not observed; closed-form unmatched branch |
| `binegbin()` / `binegbin_stanvars()` | Joint bivariate negative-binomial (overdispersed margins), a dispersion per private component |
| `binegbin_partialobs()` / `binegbin_partialobs_stanvars()` | The same, admitting rows where $y_1$ was not observed |
| `binegbin_mfd_to_dpars()` / `binegbin_dpars_to_mfd()` | Convert between $(M, f, \delta)$ and the rate/dispersion dpars |

The two constructors of a family return the same `custom_family` name, so they
share one likelihood and one set of interface functions; `*_partialobs_stanvars()`
returns exactly the same object as the plain one and exists only so that a call
reads consistently.

The `log_lik_`, `posterior_predict_` and `posterior_epred_` functions are
located by `brms` via name convention and are not normally called directly.

## References

Bürkner P-C (2017) brms: an R package for Bayesian multilevel models using Stan.
*Journal of Statistical Software* 80:1–28.

Genest C, Nešlehová J (2007) A primer on copulas for count data.
*ASTIN Bulletin* 37:475–515.

Holgate P (1964) Estimation for the Bivariate Poisson Distribution.
*Biometrika* 51:241–245. (The trivariate-reduction construction underlying these
families.)

Karlis D, Ntzoufras I (2003) Analysis of Sports Data by Using Bivariate
Poisson Models. *Journal of the Royal Statistical Society: Series D (The
Statistician)* 52:381–393.

Karlis D, Ntzoufras I (2006) Bayesian Analysis of the Differences of Count
Data. *Statistics in Medicine* 25:1885–1905.

Karlis D, Michels R, Ötting M (2026) Modelling Handball Outcomes Using
Univariate and Bivariate Approaches. *Statistical Methods & Applications*
35:263–284.

Skellam JG (1946) The Frequency Distribution of the Difference Between Two
Poisson Variates Belonging to Different Populations. *Journal of the Royal
Statistical Society* 109:296.

Smith ANH, Acuña-Marrero D, Salinas-de-León P, Harvey ES, Pawley MDM,
Anderson MJ (2020) Instantaneous versus non-instantaneous diver-operated
stereo-video (DOV) surveys of highly mobile sharks in the Galápagos Marine
Reserve. *Marine Ecology Progress Series* 649:111–123. (Supplement 3 gives the
lognormal prior grid referred to above.)

Simpson D, Rue H, Riebler A, Martins TG, Sørbye SH (2017) Penalising model
component complexity: a principled, practical approach to constructing priors.
*Statistical Science* 32:1–28.
