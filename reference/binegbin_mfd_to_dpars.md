# Convert (M, f, delta) coordinates to native binegbin/bipois dpars

Maps the interpretable coordinates – the midpoint `M` of the two
expected counts, congruence `f`, and source bias `delta` – onto the rate
dpars every family in this package takes (`mu`, `lambdaone`,
`lambdatwo`), optionally converting overdispersions to the
`shapes`/`shapexone`/`shapextwo` dpars.

The three rates are common to
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
and
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md),
so this direction serves all four constructors. The four constructors
differ in their dispersions: the Poisson families have none, and both
negative-binomial constructors have one per margin.

Everything this returns is named for a dpar a shipping family accepts,
so the output can go straight into a
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call or a
simulation. That is why `kappax` writes `shapexone` and `shapextwo` at a
common value rather than a single `shapex`: no family has had a `shapex`
dpar since 0.10.0.

[`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md)
is the inverse. It still accepts `shapex`, because its input is a stored
fit and pre-0.10.0 `binegbin` fits declare that name.

## Usage

``` r
binegbin_mfd_to_dpars(
  M,
  f,
  delta = 0,
  kappas = NULL,
  kappax = NULL,
  kappaxone = NULL,
  kappaxtwo = NULL
)
```

## Arguments

- M:

  Midpoint of the two expected counts: `mu + (lambdaone + lambdatwo)/2`,
  equal to `(E[y1] + E[y2])/2`. Non-negative.

- f:

  Congruence, the mean of the shared component as a proportion of `M`:
  `mu / M`. In `[0, 1]`. `f = 1` means perfect congruence (both
  source-exclusive components vanish); `f = 0` means no shared component
  at all.

- delta:

  Source bias on the log-ratio scale, `0.5 * log(lambdaone/lambdatwo)`.
  `0` is unbiased. `+/-Inf` is permitted and gives the limit where one
  exclusive rate is zero.

- kappas, kappax:

  Optional overdispersions, `kappa = 1/sqrt(phi)`. `0` is the Poisson
  limit. `kappas` is the overdispersion of the shared component, and the
  returned list gains `shapes` (`= 1/kappa^2`, so `kappa = 0` gives
  `Inf`). `kappax` is the shorthand for a single overdispersion
  governing *both* margins: supply it and the returned list gains
  `shapexone` and `shapextwo` at that common value, which is the
  symmetric model
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  reaches by tying the two with
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html).
  Omit both for
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
  which have no dispersion parameters.

- kappaxone, kappaxtwo:

  Optional per-margin overdispersions of the exclusive components, for
  the general case in which the two margins are free to differ. If
  supplied, the returned list gains `shapexone`/`shapextwo`. Mutually
  exclusive with `kappax`, which writes the same two slots.

## Value

A named list of `mu`, `lambdaone`, `lambdatwo`, plus `shapes` when
`kappas` is supplied, and `shapexone`/`shapextwo` when either `kappax`
or `kappaxone`/`kappaxtwo` are supplied. Every name is a dpar of a
shipping family.

## Details

Arguments are recycled to a common length, so this vectorises over
posterior draws.

**Boundary behaviour.** At `f = 1` both exclusive rates are exactly `0`
regardless of `delta` – the bias becomes unidentifiable, which
[`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md)
reports back as `NA`. This direction is always well defined; only the
inverse degenerates.

**The two spellings of the overdispersion of the exclusive components.**
Both negative-binomial constructors take the pair
`shapexone`/`shapextwo`, so this function returns `shapexone` and
`shapextwo` and never `shapex`. Supply `kappaxone`/`kappaxtwo` to give
the two margins different values, or `kappax` to give both margins a
single value – the symmetric model, term for term the pre-0.8.0
likelihood, which a fit reaches by routing both dpars through one
non-linear parameter. The two spellings are mutually exclusive in a
single call because they write the same two slots.

Before 0.10.0, `kappax` returned a dpar named `shapex` and
`kappaxone`/`kappaxtwo` returned the pair, because two different
families wanted two different things. There is now one family and one
dpar set, so both spellings produce it.

**These coordinates under the Poisson families.**
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
take the same three rates and no dispersion, so call this with `M`, `f`
and `delta` alone and pass the result straight through. There is no
`kappa` to supply: the Poisson case is not an dispersion set to a
particular value but the absence of the parameter, which is precisely
why fitting it wants its own family rather than
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
with `kappa` driven to `0`.

## See also

[`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md),
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)

## Examples

``` r
# A moderately congruent pair, source 1 recording more of the unshared
# events
binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 4.741606
#> 
#> $lambdatwo
#> [1] 3.178394
#> 

# Perfect congruence: both exclusive components vanish
binegbin_mfd_to_dpars(M = 12, f = 1, delta = 0.5)
#> $mu
#> [1] 12
#> 
#> $lambdaone
#> [1] 0
#> 
#> $lambdatwo
#> [1] 0
#> 

# kappax: one overdispersion for both exclusive components
binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4, kappax = 0.9)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 3.96
#> 
#> $lambdatwo
#> [1] 3.96
#> 
#> $shapes
#> [1] 6.25
#> 
#> $shapexone
#> [1] 1.234568
#> 
#> $shapextwo
#> [1] 1.234568
#> 

# kappaxone, kappaxtwo: one overdispersion per exclusive component
binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4,
                      kappaxone = 0.9, kappaxtwo = 0.3)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 3.96
#> 
#> $lambdatwo
#> [1] 3.96
#> 
#> $shapes
#> [1] 6.25
#> 
#> $shapexone
#> [1] 1.234568
#> 
#> $shapextwo
#> [1] 11.11111
#> 

# To FIT in these coordinates, pass them through a non-linear formula
# (every dpar is log-linked, so the link supplies the exp()):
#   bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#     nlf(mu        ~ logM + log_inv_logit(logitf)) +
#     nlf(lambdaone ~ log(2) + logM + log_inv_logit(-logitf) +
#                     log_inv_logit(2 * delta)) +
#     nlf(lambdatwo ~ log(2) + logM + log_inv_logit(-logitf) +
#                     log_inv_logit(-2 * delta)) +
#     lf(logM ~ 1, logitf ~ 1, delta ~ 1)
# where logM = log M, logitf = logit f and delta = artanh(beta).
```
