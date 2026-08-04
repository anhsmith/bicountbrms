# Deprecated names for the censoring-aware Negative-Binomial family

\`binegbin_joint()\` and \`binegbin_joint_stanvars()\` were renamed to
\[binegbin_cens()\] and \[binegbin_cens_stanvars()\] in 0.9.0. The old
names still work and return the new objects, with a deprecation warning.

\`log_lik_binegbin_joint()\`, \`posterior_predict_binegbin_joint()\` and
\`posterior_epred_binegbin_joint()\` exist so that fits made before
0.9.0 keep working. Such a fit stores \`family\$name ==
"binegbin_joint"\`, and brms builds its post-processing method names
from that stored name, so these three are reached by \`loo()\`,
\`posterior_predict()\` and \`posterior_epred()\` without the user
naming them. They forward silently – the name is a property of the
stored fit, not of anything the caller can change, and brms calls the
first two once per observation. \*\*No refitting is required.\*\*

New code should use \[binegbin_cens()\]. These names will be removed in
the next major version.

## Usage

``` r
binegbin_joint()

binegbin_joint_stanvars()

log_lik_binegbin_joint(i, prep)

posterior_predict_binegbin_joint(i, prep, ...)

posterior_epred_binegbin_joint(prep)
```

## Arguments

- i, prep:

  Passed through unchanged to the corresponding \`binegbin_cens\`
  method.

- ...:

  Passed through unchanged.

## Value

As the corresponding \[binegbin_cens()\] function.
