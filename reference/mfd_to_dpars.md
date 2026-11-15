# Convert from the congruence–bias parameterisation to native parameters

Converts from the congruence–bias parameterisation (the midpoint `M` of
the two expected counts, the congruence `f`, and the source bias
`delta`) to the native parameters (the rates `mu`, `lambdaone`, and
`lambdatwo`), and, optionally, converts overdispersions (`kappas`,
`kappax`, `kappaxone`, and `kappaxtwo`) to the native dispersions
(`shapes`, `shapexone`, and `shapextwo`).
[`dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/dpars_to_mfd.md)
converts in the other direction.

## Usage

``` r
mfd_to_dpars(
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
  exclusive rates are zero); `f = 0` means no shared component.

- delta:

  Source bias on the log-ratio scale, `0.5 * log(lambdaone/lambdatwo)`.
  `0` is unbiased. `+/-Inf` is permitted and gives the limit where one
  exclusive rate is zero.

- kappas, kappax:

  Optional overdispersions on the scale `kappa = 1/sqrt(phi)`, where
  `phi` is the NB2 (inverse) dispersion parameter; `kappa = 0` is the
  Poisson limit. `kappas` is the overdispersion of the shared component.
  Supplying `kappas` adds `shapes = 1/kappas^2` to the returned list, so
  `kappas = 0` gives `shapes = Inf`. `kappax` is one overdispersion for
  both exclusive components. Supplying `kappax` adds `shapexone` and
  `shapextwo` to the returned list, both equal to `1/kappax^2`. For
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
  which have no dispersion parameters, leave every overdispersion
  argument unset.

- kappaxone, kappaxtwo:

  Optional overdispersions of the two exclusive components, for the
  general case in which they differ. Supplying them adds `shapexone` and
  `shapextwo` to the returned list. Mutually exclusive with `kappax`.

## Value

A named list of `mu`, `lambdaone`, and `lambdatwo`, with `shapes` added
when `kappas` is supplied, and `shapexone` and `shapextwo` added when
`kappax`, or `kappaxone` and `kappaxtwo`, are supplied.

## Details

Each argument may be a vector, such as the posterior draws of `M`, `f`,
and `delta`. `M`, `f`, and `delta` are recycled to the length of the
longest of the three, and any overdispersions to the same length; each
element of the returned list is a vector of that length.

To fit a model in the congruence–bias parameterisation, use
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md),
as in
[`vignette("bicountbrms")`](https://anhsmith.github.io/bicountbrms/articles/bicountbrms.md).

## See also

[`dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/dpars_to_mfd.md),
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)

## Examples

``` r
# A moderately congruent pair, source 1 recording more of the exclusive
# events
mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 4.741606
#> 
#> $lambdatwo
#> [1] 3.178394
#> 

# Perfect congruence: both exclusive rates are zero
mfd_to_dpars(M = 12, f = 1, delta = 0.5)
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
mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4, kappax = 0.9)
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
mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4,
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
```
