# Discrete-Laplace custom family for brms (location 0, free scale)

Returns a brms custom family for the discrete Laplace distribution,
location fixed at 0, discretised from the continuous Laplace(0, b) via
CDF differencing: \`P(Z=z) = F(z+0.5) - F(z-0.5)\`. One parameter, sigma
(link = "log"), the SD; the mean is always zero. Unlike skellam1/
skellam2, the PMF and CCDF are closed-form
(\`double_exponential_lcdf\`-based – Stan's name for the Laplace
distribution – no Bessel function, no large-argument branch or iteration
cap needed).

Use in a brm() call as: brm(y ~ ..., family = dlaplace1(), stanvars =
dlaplace1_stanvars(), data = ...)

## Usage

``` r
dlaplace1()

dlaplace1_stanvars()

log_lik_dlaplace1(i, prep)

posterior_predict_dlaplace1(i, prep, ...)

posterior_epred_dlaplace1(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Naming note.\*\* Same forced naming as \`skellam1()\`:
\`brms::custom_family()\` requires a dpar literally named \`"mu"\`; here
it represents sigma (the SD), not a mean. See \`?skellam1\` Details for
the full rationale.

\*\*sigma-to-b conversion.\*\* Stan's \`double_exponential_lcdf\`
expects the continuous Laplace's own scale parameter, \`b\`.
Var(Laplace(0,b)) = \`2\*b^2\`, so SD = \`b\*sqrt(2)\`; treating sigma
as exactly that SD (the discretisation perturbs the true discrete
variance only slightly, and this keeps sigma on the same scale as the
other three families) gives \`b = sigma / sqrt(2)\`, computed first in
both \`dlaplace1_lpmf\` and \`dlaplace1_lccdf\`.

\*\*Validation note.\*\* \`extraDistr::ddlaplace()\` implements a
different discrete Laplace — its \`scale\` argument is actually a decay
probability \`p\` for the exact closed form \`P(z) = (1-p)/(1+p) \*
p^\|z\|\`, not a continuous-Laplace \`b\` — confirmed numerically to NOT
match this family's CDF-differenced PMF (e.g. at \`b=3\`,
\`p=exp(-1/3)\`: \`P(0) = 0.1535\` here vs \`0.1651\` there). This
package's tests validate against a hand-derived CDF-difference R
reference instead — the documented fallback for when a package reference
isn't applicable.
