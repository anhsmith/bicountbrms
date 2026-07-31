# Symmetric Skellam custom family for brms

Returns a brms custom family for the symmetric Skellam distribution,
Skellam(mu_skellam, mu_skellam) — the distribution of the difference of
two independent Poisson(mu_skellam) random variables. The single
parameter is sigma (link = "log"), the SD of that difference; the mean
is always zero. Internally, mu_skellam = sigma^2 / 2 is derived as a
transformed parameter and fed to the underlying Bessel-function PMF,
which is otherwise unchanged.

Use in a brm() call as: brm(y ~ ..., family = skellam1(), stanvars =
skellam1_stanvars(), data = ...)

## Usage

``` r
skellam1()

skellam1_stanvars()

log_lik_skellam1(i, prep)

posterior_predict_skellam1(i, prep, ...)

posterior_epred_skellam1(prep)
```

## Value

A brms custom_family object.

## Details

This family was originally parameterised directly on mu_skellam (link =
"log"); it now samples on sigma instead, for a common (mean, SD-scale)
convention shared with skellam2(), dlaplace1(), and dlaplace2(). Since
sigma = sqrt(2 \* mu_skellam), a prior previously stated on
log(mu_skellam) — e.g. normal(1, 1.5) — translates as: log(sigma) = 0.5
\* log(2) + 0.5 \* log(mu_skellam) so an intercept of 1 on the old
log(mu_skellam) scale corresponds to an intercept of 0.5\*log(2) +
0.5\*1 ≈ 0.847 on the new log(sigma) scale, and the old prior's SD of
1.5 becomes 0.75 on the new scale (a linear transform of a normal is
normal). This is a scale correspondence only — slope-coefficient
interpretations from the old parameterisation are NOT carried forward;
any offset-vs-free-slope diagnostic should be redone fresh against this
sigma-scale parameterisation.

\*\*Naming note.\*\* \`brms::custom_family()\` hard-requires one
\`dpars\` entry to be literally named \`"mu"\` (\`stop2("All families
must have a 'mu' parameter.")\`, unconditional, no override) — every
family built on it, including this one, must comply regardless of what
that parameter actually represents. For skellam1, the brms/Stan-level
dpar named \`mu\` IS sigma (the SD of the difference, log-linked); it is
NOT the distribution's mean, which is structurally zero throughout. This
is a forced naming collision with brms's API, not a reversion to the
pre-reparameterisation behaviour: internally, \`mu_skellam = mu^2 / 2\`
is still derived from it before reaching the Bessel-function PMF,
exactly as documented above for "sigma". All R-side helper functions
below immediately rebind this dpar to a variable called \`sigma\` so
that no code in this package, other than the literal \`dpars\`/
\`get_dpar()\` calls forced by brms, ever refers to it as \`mu\`.
