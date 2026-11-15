# Convert from native parameters to the congruence–bias parameterisation

Converts from the native parameters (the rates `mu`, `lambdaone`, and
`lambdatwo`) to the congruence–bias parameterisation (the midpoint `M`
of the two expected counts, the congruence `f`, and the source bias
`delta`), and, optionally, converts native dispersions (`shapes`,
`shapex`, `shapexone`, and `shapextwo`) to overdispersions (`kappas`,
`kappax`, `kappaxone`, and `kappaxtwo`).
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
converts in the other direction.

## Usage

``` r
dpars_to_mfd(
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

  Rate of the shared component.

- lambdaone, lambdatwo:

  Rates of the source-1-exclusive and source-2-exclusive components.

- shapes, shapex:

  Optional NB2 (inverse) dispersion parameters. Supplying `shapes` adds
  `kappas = 1/sqrt(shapes)` to the returned list, so `shapes = Inf`
  gives `kappas = 0`. A `binegbin` fit made before 0.10.0 has one
  dispersion, `shapex`, for both exclusive components; supplying
  `shapex` adds `kappax`. For
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
  which have no dispersion parameters, leave every dispersion argument
  unset.

- shapexone, shapextwo:

  Optional NB2 (inverse) dispersion parameters of the two exclusive
  components. Supplying them adds `kappaxone` and `kappaxtwo` to the
  returned list. Mutually exclusive with `shapex`.

## Value

A named list of `M`, `f`, `delta`, and `beta = tanh(delta)`, the source
bias on a scale bounded by -1 and 1, with `kappas`, `kappax`,
`kappaxone`, and `kappaxtwo` added when the corresponding dispersions
are supplied.

## Details

Each argument may be a vector, such as the posterior draws of the rates.
`mu`, `lambdaone`, and `lambdatwo` are recycled to the length of the
longest of the three, and any dispersions to the same length; each
element of the returned list is a vector of that length.

Rates of exactly zero, which posterior draws do not produce, give values
that are not finite: `delta` is `NA` when both exclusive rates are zero,
and `Inf` or `-Inf` when one of them is, and `f` is `NA` when `M` is
zero.

## See also

[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md),
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md),
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)

## Examples

``` r
dpars_to_mfd(mu = 8.04, lambdaone = 4.75, lambdatwo = 3.17)
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
d <- mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)[c("M", "f", "delta")]
#> $M
#> [1] 12
#> 
#> $f
#> [1] 0.67
#> 
#> $delta
#> [1] 0.2
#> 

# Perfect congruence: delta cannot be recovered and is returned as NA
dpars_to_mfd(mu = 12, lambdaone = 0, lambdatwo = 0)$delta
#> [1] NA
```
