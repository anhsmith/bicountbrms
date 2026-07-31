# Joint bivariate-Poisson custom family for brms

Returns a brms custom family for the joint distribution of a matched
pair of counts, \`(y1, y2)\`, constructed via trivariate reduction: \`y1
= N_shared + N1\`, \`y2 = N_shared + N2\`, with \`N_shared ~
Poisson(mu)\`, \`N1 ~ Poisson(lambdaone)\`, \`N2 ~ Poisson(lambdatwo)\`
mutually independent given their rates. All three rates are link =
"log". Modelling the pair jointly avoids regressing the difference on
one of its own components (a \`d = y1 - y2 ~ y2\` design), which induces
regression to the mean.

\`y1\` is the family's response; \`y2\` is passed in as supplementary
integer data via brms's \`vint()\` addition term, since brms's
\`custom_family()\` machinery is built around a single declared response
column – see Details.

Use in a brm() call as: brm( bf(y1 \| vint(y2) ~ ...), family =
bipois(), stanvars = bipois_stanvars(), data = dat )

## Usage

``` r
bipois()

bipois_stanvars()

log_lik_bipois(i, prep)

posterior_predict_bipois(i, prep, ...)

posterior_epred_bipois(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Naming note.\*\* Same forced naming as
\`skellam1()\`/\`dlaplace1()\`/ \`dnorm1()\`: \`brms::custom_family()\`
requires a dpar literally named \`"mu"\` (\`stop2("All families must
have a 'mu' parameter.")\`, unconditional). Here it is bound to
\`lambda_shared\`, the rate of the component shared between \`y1\` and
\`y2\` – not a mean of either response individually. \`lambdaone\`
(source-1-only rate) and \`lambdatwo\` (source-2-only rate) are the
other two dpars, plainly named (no forced reinterpretation needed for
those two).

\*\*Why \`lambdaone\`, not \`lambda1\`.\*\* \`custom_family()\` rejects
dpar names ending in a digit (\`stop2("'dpars' should not end with a
number.")\`), as well as dots and underscores. The documentation
therefore writes these rates as \\\lambda_1\\ and \\\lambda_2\\ while
the code must spell them \`lambdaone\`/\`lambdatwo\`. See the notation
table in the package README.

\*\*Why \`y2\` travels via \`vint()\`, not as a second response.\*\*
brms's \`custom_family()\` API supports exactly one declared response
column (\`Y\`) plus optional supplementary integer/real data
(\`vint()\`/\`vreal()\` addition terms) – the same mechanism used for,
e.g., binomial trial counts in the brms custom-families vignette. There
is no \*undeclared-response\* concept for a genuinely joint two-count
likelihood; \`vint(y2)\` is the correct fit for that gap, not a
workaround. This does mean \`y2\` is \*not\* itself treated as
brms-modelled response data (no missing-value handling, no resp\_\*()
addition terms apply to it) – it is fixed, observed per-row data,
consistent with the fact that every row used here comes from the matched
(both-observed) subset.

\*\*Order of dpars matters for the generated Stan call.\*\* brms
generates \`target += bipois_lpmf(Y\[n\] \| mu\[n\], lambdaone\[n\],
lambdatwo\[n\], vint1\[n\])\` – dpars in the order declared here, then
vint/vreal args in the order declared in \`vars\`. \`bipois_stan_funs\`
(stanfunctions via \`bipois_stanvars()\`) declares \`bipois_lpmf\` with
exactly this argument order; changing the order here without changing
the Stan signature (or vice versa) silently swaps which rate governs
which count.
