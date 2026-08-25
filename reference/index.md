# Package index

## Joint bivariate-count families

Model the matched pair jointly via trivariate reduction, capturing its
correlation, marginal overdispersion, and difference together rather
than the difference alone. Each component distribution has one family
name and two constructors: the plain one for a fully paired design, and
the \_partialobs one for a design in which the first count is missing on
some rows. Both share a likelihood and a set of post-processing methods.

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
  : Joint bivariate-Negative-Binomial custom family for brms

## Partially observed pairs

The same two models, for rows on which the first count was never
recorded. Those rows are scored by the integrated-out marginal of the
same joint, so they still inform the shared component and the second
source rather than being dropped. For Poisson components that marginal
is closed form; for Negative-Binomial components it is a convolution.
This is unrelated to brms’s own cens() addition term, which means a
value known to lie in a set.

- [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  [`bipois_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
  : Joint bivariate-Poisson family for partially observed pairs
- [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  [`binegbin_partialobs_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
  : Joint bivariate-Negative-Binomial family for partially observed
  pairs

## Parameterisation helpers

Convert between the dpars a family takes and the interpretable
coordinates of overall level, congruence, and source bias.

- [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md)
  : Convert (M, f, delta) coordinates to native binegbin/bipois dpars
- [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md)
  : Convert native binegbin/bipois dpars to (M, f, delta) coordinates
