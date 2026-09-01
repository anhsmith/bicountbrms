# Convert native binegbin/bipois dpars to (M, f, delta) coordinates

Inverse of
[`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md).
Reads the rate dpars `mu`, `lambdaone`, `lambdatwo` back into the
interpretable overall level `M`, congruence `f`, and source bias
`delta`, optionally converting NB2 dispersions back to the SD scale.

As with the forward direction, the three rates are common to
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
and
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md),
so this serves all four constructors; only the dispersion arguments
differ.

This direction reads a *stored fit*, so it accepts `shapex` – the single
excess dispersion declared by every `binegbin` fit made before 0.10.0 –
as well as the current `shapexone`/`shapextwo`. The forward direction
writes only current names; see its documentation for why the two differ.

## Usage

``` r
binegbin_dpars_to_mfd(
  mu,
  lambdaone,
  lambdatwo,
  shapes = NULL,
  shapex = NULL,
  shapexone = NULL,
  shapextwo = NULL
)
```

## Arguments

- mu:

  Shared-component rate.

- lambdaone, lambdatwo:

  The two excess rates.

- shapes, shapex:

  Optional NB2 dispersions. If supplied, the returned list gains
  `kappas`/`kappax` (`= 1/sqrt(shape)`, so `shape = Inf` gives `0`).
  `shapex` is the single excess dispersion declared by a `binegbin` fit
  made before 0.10.0, when one dpar governed both margins; no shipping
  family takes it now. A fit of
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  or
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  has no dispersion to pass: read its three rates alone.

- shapexone, shapextwo:

  Optional per-margin NB2 excess dispersions, as taken by both
  negative-binomial constructors since 0.10.0. If supplied, the returned
  list gains `kappaxone`/`kappaxtwo`. Mutually exclusive with `shapex`,
  which is the same quantity under its older name.

## Value

A named list of `M`, `f`, `delta`, plus `beta` (the bounded bias
`tanh(delta)`), `kappas`/`kappax` when `shapes`/`shapex` are supplied,
and `kappaxone`/`kappaxtwo` when `shapexone`/`shapextwo` are.

## Details

Arguments are recycled to a common length, so this vectorises over
posterior draws – e.g. to convert a whole posterior into interpretable
coordinates.

**Boundary behaviour**, which the forward direction does not have:

- `lambdaone == lambdatwo == 0` (perfect congruence, `f = 1`): the bias
  is genuinely unidentified – there is no excess to be biased – and
  `delta` is returned as `NA`, not `0`. Zero would assert an unbiased
  source, which the data at that point cannot support.

- `M == 0` (nothing anywhere): `f` is undefined and returned as `NA`.

- Exactly one excess rate `0`: `delta` is `+/-Inf`, the well-defined
  limit where one source never records an unshared event.

Because of the first case, round-tripping is exact everywhere except at
`f = 1`, where `delta` cannot be recovered. Hold one coordinate system
as the source of truth rather than repeatedly converting back and forth.

## See also

[`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md),
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)

## Examples

``` r
binegbin_dpars_to_mfd(mu = 8.04, lambdaone = 4.75, lambdatwo = 3.17)
#> $M
#> [1] 12
#> 
#> $f
#> [1] 0.67
#> 
#> $delta
#> [1] 0.2022065
#> 
#> $beta
#> [1] 0.1994949
#> 

# Round trip
d <- binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)[c("M", "f", "delta")]
#> $M
#> [1] 12
#> 
#> $f
#> [1] 0.67
#> 
#> $delta
#> [1] 0.2
#> 

# Perfect congruence: bias is unidentified, reported as NA
binegbin_dpars_to_mfd(mu = 12, lambdaone = 0, lambdatwo = 0)$delta
#> [1] NA
```
