# Package index

## Joint bivariate-count families

Model the matched pair jointly via trivariate reduction, capturing its
correlation, marginal overdispersion, and difference together rather
than the difference alone.

- [`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  [`bipois_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  [`log_lik_bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  [`posterior_predict_bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  [`posterior_epred_bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md)
  : Joint bivariate-Poisson custom family for brms
- [`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  [`binegbin_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  [`log_lik_binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  [`posterior_predict_binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  [`posterior_epred_binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  : Joint bivariate-Negative-Binomial custom family for brms
- [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
  [`binegbin_joint_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
  [`log_lik_binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
  [`posterior_predict_binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)
  : Censoring-aware joint bivariate-Negative-Binomial family for brms

## Difference families

Model the difference `d = y1 - y2` directly. The `1` variants fix the
location at zero (do the two sources agree on average?); the `2`
variants estimate it (by how much do they disagree?).

- [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  [`skellam1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  [`log_lik_skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  [`posterior_predict_skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  [`posterior_epred_skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)
  : Symmetric Skellam custom family for brms
- [`skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
  [`skellam2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
  [`log_lik_skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
  [`posterior_predict_skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
  [`posterior_epred_skellam2()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2.md)
  : Asymmetric Skellam custom family for brms
- [`dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  [`dnorm1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  [`log_lik_dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  [`posterior_predict_dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  [`posterior_epred_dnorm1()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1.md)
  : Discrete-normal custom family for brms (location 0, free scale)
- [`dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
  [`dnorm2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
  [`log_lik_dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
  [`posterior_predict_dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
  [`posterior_epred_dnorm2()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2.md)
  : Discrete-normal custom family for brms (free location and scale)
- [`dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  [`dlaplace1_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  [`log_lik_dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  [`posterior_predict_dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  [`posterior_epred_dlaplace1()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1.md)
  : Discrete-Laplace custom family for brms (location 0, free scale)
- [`dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  [`dlaplace2_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  [`log_lik_dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  [`posterior_predict_dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  [`posterior_epred_dlaplace2()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2.md)
  : Discrete-Laplace custom family for brms (free location and scale)

## Truncation support

Stan `lccdf` functions injected as stanvars so the difference families
normalise correctly under
[`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html).

- [`skellam1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1_lccdf_stanvars.md)
  : Truncated-Skellam log-CCDF for use with brms's resp_trunc()
- [`skellam2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2_lccdf_stanvars.md)
  : Truncated-asymmetric-Skellam log-CCDF for use with brms's
  resp_trunc()
- [`dnorm1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm1_lccdf_stanvars.md)
  : Truncated-discrete-normal log-CCDF for use with brms's resp_trunc()
- [`dnorm2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dnorm2_lccdf_stanvars.md)
  : Truncated-discrete-normal log-CCDF for use with brms's resp_trunc()
  (free location and scale)
- [`dlaplace1_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace1_lccdf_stanvars.md)
  : Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()
- [`dlaplace2_lccdf_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/dlaplace2_lccdf_stanvars.md)
  : Truncated-discrete-Laplace log-CCDF for use with brms's resp_trunc()
  (free location and scale)

## Parameterisation helpers

Convert between the dpars a family takes and more interpretable
coordinates.

- [`skellam2_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam2_dpars.md)
  : Report skellam2's derived quantities from a fitted model
- [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  : Convert (M, f, delta) coordinates to native binegbin/bipois dpars
- [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_dpars_to_mfd.md)
  : Convert native binegbin/bipois dpars to (M, f, delta) coordinates
