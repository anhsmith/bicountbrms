# Package index

## Joint bivariate-count families

Model the matched pair jointly via trivariate reduction, capturing its
correlation, marginal overdispersion, and difference together rather
than the difference alone. Both counts must be observed on every row.

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

## Censoring-aware families

The same two models, extended to rows on which the first count is
unobserved, through the integrated-out marginal of the same joint. For
Poisson components that marginal is closed form; for Negative-Binomial
components it is a convolution, and the two source-specific components
may differ in overdispersion.

- [`bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  [`bipois_cens_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  [`log_lik_bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  [`posterior_predict_bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  [`posterior_epred_bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  : Censoring-aware joint bivariate-Poisson family for brms
- [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  [`binegbin_cens_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  [`log_lik_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  [`posterior_predict_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  [`posterior_epred_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  : Censoring-aware joint bivariate-Negative-Binomial family for brms

## Parameterisation helpers

Convert between the dpars a family takes and the interpretable
coordinates of overall level, congruence, and source bias.

- [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md)
  : Convert (M, f, delta) coordinates to native binegbin/bipois dpars
- [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md)
  : Convert native binegbin/bipois dpars to (M, f, delta) coordinates

## Deprecated

Retained so that code and fitted models predating the 0.9.0 rename
binegbin_joint -\> binegbin_cens keep working. Removed in the next major
version.

- [`binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  [`binegbin_joint_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  [`log_lik_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  [`posterior_predict_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  [`posterior_epred_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  : Deprecated names for the censoring-aware Negative-Binomial family
