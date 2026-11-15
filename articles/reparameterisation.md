# The congruence–bias parameterisation

The widget below demonstrates the recommended congruence–bias
parameterisation by simulating ten pairs of counts from a model (see
[`rbicount()`](https://anhsmith.github.io/bicountbrms/reference/rbicount.md)).
The relationship with the native parameterisation is shown by moving the
sliders: a change to a parameter in one parameterisation changes the
parameters in the other. The plot below the sliders shows the expected
decomposition of the two counts (left) and the simulated pairs (right).

re-simulate

seed

reset

#### Congruence–bias parameterisation

#### Native parameterisation

The coordinates of the congruence–bias parameterisation directly
represent the quantities usually of interest in a comparison of two
sources:

| Coordinate | Meaning | Definition | Name in brms |
|----|----|----|----|
| M | midpoint of the expected counts of the two sources | M = \mu + (\lambda_1 + \lambda_2)/2 | `logM`, for \log M |
| f | congruence, the proportion of M recorded by both sources | f = \mu / M | `logitf`, for \operatorname{logit} f |
| \beta | source bias, zero when the two exclusive rates are equal | \beta = (\lambda_1 - \lambda_2)/(\lambda_1 + \lambda_2) | `delta`, for \operatorname{artanh}\beta |
| \kappa_s, \kappa\_{x1}, \kappa\_{x2} | overdispersions | \kappa = 1/\sqrt{\phi} | `kappas`, `kappaxone`, `kappaxtwo` |

The native parameters are the rates and dispersions of the three latent
components:

| Parameter | Meaning | Definition | Name in brms |
|----|----|----|----|
| \mu | rate of the shared component, recorded by both sources | \mu = Mf | `mu` |
| \lambda_1 | rate of the source-exclusive component of source 1 | \lambda_1 = M(1 - f)(1 + \beta) | `lambdaone` |
| \lambda_2 | rate of the source-exclusive component of source 2 | \lambda_2 = M(1 - f)(1 - \beta) | `lambdatwo` |
| \phi_s, \phi\_{x1}, \phi\_{x2} | NB2 (inverse) dispersions | \phi = 1/\kappa^2 | `shapes`, `shapexone`, `shapextwo` |

The functions
[`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
and
[`dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/dpars_to_mfd.md)
convert between the two parameterisations, for single values or for
vectors of posterior draws:

``` r

library(bicountbrms)

dpars <- mfd_to_dpars(
  M      = 12,
  f      = 0.67,
  delta  = atanh(0.2),
  kappas = 0.6,
  kappax = 1.0
)
unlist(dpars)
#>        mu lambdaone lambdatwo    shapes shapexone shapextwo 
#>  8.040000  4.752000  3.168000  2.777778  1.000000  1.000000

mfd <- dpars_to_mfd(
  mu        = dpars$mu,
  lambdaone = dpars$lambdaone,
  lambdatwo = dpars$lambdatwo
)
unlist(mfd[c("M", "f", "beta")])
#>     M     f  beta 
#> 12.00  0.67  0.20
```
