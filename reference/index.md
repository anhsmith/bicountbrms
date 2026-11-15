# Package index

## Families for paired counts

Constructors for the bipois and binegbin families, which model a pair of
counts by trivariate reduction with Poisson or negative-binomial
components. Each family has a second constructor for partially paired
designs, listed below, and the two constructors of a family share a
likelihood and the post-processing methods.

- [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  [`bipois_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  [`log_lik_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  [`posterior_predict_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  [`posterior_epred_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  : Joint bivariate-Poisson custom family for brms
- [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  [`binegbin_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  [`log_lik_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  [`posterior_predict_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  [`posterior_epred_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  : Joint bivariate custom family with negative-binomial components for
  brms

## Partially paired designs

Constructors for partially paired designs, in which source 1 has missing
values. On an unpaired row, the likelihood term is the marginal
probability of y2, so the row informs the shared component and the
exclusive component of source 2.

- [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  [`bipois_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  : Joint bivariate-Poisson family for partially observed pairs
- [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  [`binegbin_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  : Joint bivariate family with negative-binomial components for
  partially observed pairs

## Prediction and imputation

Draws and expectations of both counts together, and of y1 given the
recorded y2. The functions conditional on y2 impute y1 on the unpaired
rows of a fit from bipois_partialobs() or binegbin_partialobs().

- [`posterior_predict_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
  [`posterior_epred_y1_and_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_and_y2.md)
  : Draws and expectations of both counts from the joint model
- [`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
  [`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
  : Draws and expectations of y1 given y2

## Congruence–bias coordinates

Conversion between the native parameters and the congruence–bias
coordinates (the midpoint of the two expected counts, the congruence,
and the source bias), and the brms formula for a model in the
coordinates.

- [`mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/mfd_to_dpars.md)
  : Convert from the congruence–bias parameterisation to native
  parameters
- [`dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/dpars_to_mfd.md)
  : Convert from native parameters to the congruence–bias
  parameterisation
- [`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md)
  : Build a brms formula in the congruence–bias parameterisation

## Simulation

Random generation of paired counts from the generating model of the
bipois and binegbin families, with the shared and source-exclusive
components of each pair.

- [`rbicount()`](https://anhsmith.github.io/bicountbrms/reference/rbicount.md)
  : Random generation of paired counts
