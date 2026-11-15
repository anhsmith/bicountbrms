# Random generation of paired counts

`rbicount()` draws pairs of counts `(y1, y2)` from the generating model
of the `bipois` and `binegbin` families, and returns the three latent
components of each pair with the pair. The arguments are the native
parameters of the families, in the order and under the names that
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
returns them.

## Usage

``` r
rbicount(
  n,
  mu,
  lambdaone,
  lambdatwo,
  shapes = NULL,
  shapexone = NULL,
  shapextwo = NULL
)
```

## Arguments

- n:

  Number of pairs to draw. If `length(n) > 1`, the length of `n` is
  taken to be the number required, as in
  [`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html).

- mu:

  Mean of the shared component `Ns`.

- lambdaone, lambdatwo:

  Means of the source-exclusive components `N1` and `N2`.

- shapes, shapexone, shapextwo:

  NB2 (inverse) dispersion parameters of `Ns`, `N1`, and `N2`. Supply
  all three for negative-binomial components (the `binegbin` family), or
  none for Poisson components (the `bipois` family).

## Value

A data frame with `n` rows and five integer columns: `y1` and `y2`, the
two recorded counts, followed by `Ns`, `N1`, and `N2`, the shared
component and the two exclusive components, with `y1 = Ns + N1` and
`y2 = Ns + N2` on every row.

## Details

Both families construct the pair by trivariate reduction. Three
independent latent counts are drawn: the shared component `Ns`, recorded
by both sources, with mean `mu`; and the source-exclusive components
`N1` and `N2`, each recorded by one source alone, with means `lambdaone`
and `lambdatwo`. The recorded counts are

\$\$y_1 = N_s + N_1, \qquad y_2 = N_s + N_2,\$\$

so that \\E\[y_1\] = \mu + \lambda_1\\, \\E\[y_2\] = \mu + \lambda_2\\,
and \\\mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_s)\\.

With the three dispersions supplied, each component is negative binomial
in the NB2 parameterisation, drawn by `rnbinom(size = phi, mu = m)`. A
component with mean \\m\\ and dispersion \\\phi\\ has variance \\m +
m^2/\phi\\, so a small \\\phi\\ gives a heavy tail and \\\phi = \infty\\
gives the Poisson limit. The dispersion of `Ns` is `shapes`, of `N1` is
`shapexone`, and of `N2` is `shapextwo`. With all three dispersions
`NULL`, each component is Poisson, with variance equal to its mean.
Supplying some of the dispersions and not others stops with an error,
because neither family has a mixture of Poisson and negative-binomial
components.

Each parameter is recycled to length `n`, as by
[`stats::rnbinom()`](https://rdrr.io/r/stats/NegBinomial.html), so a
parameter may be a single value or a vector with one value per pair. The
list returned by
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
has elements named `mu`, `lambdaone`, `lambdatwo`, and, when
overdispersions are supplied to it, `shapes`, `shapexone`, and
`shapextwo`, so `do.call(rbicount, c(list(n = n), pars))` draws pairs at
the native parameters `pars` (see Examples). The overdispersion \\\kappa
= 0\\ gives `shapes = Inf` in
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md),
and `rbicount()` then draws a Poisson component.

The components are drawn by the code that
[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
uses for each posterior draw, so a pair simulated by `rbicount()` and a
posterior predictive draw at the same parameter values have the same
distribution.

To simulate a partially paired design, draw an observation indicator,
such as `y1_obs <- rbinom(n, 1, 0.2)`, and set `y1` to `-1` on the rows
where `y1_obs` is 0 before fitting with
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
or
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).
The columns `Ns` and `N1` keep the values from which the unrecorded `y1`
would have been formed.

## See also

[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
for the native parameters at given values of the midpoint, congruence,
and source bias;
[`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
for draws of the pair from a fitted model.

## Examples

``` r
set.seed(1)

# Poisson components, as in bipois()
rbicount(5, mu = 4, lambdaone = 1, lambdatwo = 2)
#>   y1 y2 Ns N1 N2
#> 1  5  4  3  2  1
#> 2  6  4  3  3  1
#> 3  5  7  4  1  3
#> 4  8  8  7  1  1
#> 5  2  5  2  0  3

# Negative-binomial components, as in binegbin(). The covariance of the
# two counts is Var(Ns) = 4 + 4^2/2 = 12.
d <- rbicount(
  10000,
  mu        = 4,
  lambdaone = 1,
  lambdatwo = 2,
  shapes    = 2,
  shapexone = 4,
  shapextwo = 4
)
cov(d$y1, d$y2)
#> [1] 11.58208

# Parameters from the midpoint, congruence, and source bias
pars <- mfd_to_dpars(
  M      = 5,
  f      = 0.7,
  delta  = atanh(0.3),
  kappas = 0.5,
  kappax = 0.5
)
d <- do.call(rbicount, c(list(n = 100), pars))

# A partially paired design. y2 is recorded on every unit, and y1 on about
# 20% of units.
d$y1_obs  <- rbinom(nrow(d), 1, 0.2)
d$y1_true <- d$y1
d$y1[d$y1_obs == 0] <- -1L
head(d)
#>   y1 y2 Ns N1 N2 y1_obs y1_true
#> 1 -1  6  5  0  1      0       5
#> 2 -1  4  4  1  0      0       5
#> 3 -1  8  8 10  0      0      18
#> 4 -1  6  4  3  2      0       7
#> 5 -1  3  3  0  0      0       3
#> 6 -1  4  4  1  0      0       5

# Thirty units in each of three treatments that differ in congruence
if (requireNamespace("dplyr", quietly = TRUE)) {
  units <- data.frame(
    treatment = rep(c("a", "b", "c"), each = 30)
  )
  units |>
    dplyr::mutate(
      f = c(a = 0.5, b = 0.7, c = 0.9)[treatment],
      as.data.frame(
        mfd_to_dpars(
          M      = 5,
          f      = f,
          delta  = 0,
          kappas = 0.5,
          kappax = 0.5
        )
      ),
      rbicount(
        dplyr::n(),
        mu,
        lambdaone,
        lambdatwo,
        shapes,
        shapexone,
        shapextwo
      )
    ) |>
    head()
}
#>   treatment   f  mu lambdaone lambdatwo shapes shapexone shapextwo y1 y2 Ns N1
#> 1         a 0.5 2.5       2.5       2.5      4         4         4  7  7  6  1
#> 2         a 0.5 2.5       2.5       2.5      4         4         4  5  7  3  2
#> 3         a 0.5 2.5       2.5       2.5      4         4         4  7  5  3  4
#> 4         a 0.5 2.5       2.5       2.5      4         4         4  7  2  1  6
#> 5         a 0.5 2.5       2.5       2.5      4         4         4  4  3  1  3
#> 6         a 0.5 2.5       2.5       2.5      4         4         4  5  4  3  2
#>   N2
#> 1  1
#> 2  4
#> 3  2
#> 4  1
#> 5  2
#> 6  1
```
