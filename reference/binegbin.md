# Joint bivariate-Negative-Binomial custom family for brms

Overdispersed sibling of \[bipois()\]. Returns a brms custom family for
the joint distribution of a matched count pair \`(y1, y2)\` via
trivariate reduction with Negative-Binomial (rather than Poisson) latent
components: \`y1 = N_shared + N1\`, \`y2 = N_shared + N2\`, with
\`N_shared ~ NB2(mu, shapes)\`, \`N1 ~ NB2(lambdaone, shapex)\`, \`N2 ~
NB2(lambdatwo, shapex)\` mutually independent given their rates.
\`NB2(m, phi)\` has mean \`m\` and variance \`m + m^2/phi\` (Stan
\`neg_binomial_2\`; R \`dnbinom(size = phi, mu = m)\`).

Five dpars: the three rates (\`mu\` = shared rate,
\`lambdaone\`/\`lambdatwo\` = the two source-specific rates) plus two
dispersions – \`shapes\` for the shared component and \`shapex\` shared
across the two excess components. All five use \`link = "log"\`. Supply
the excess rates through a non-linear formula without an explicit
\`exp()\` (the log link applies it): \`nlf(lambdaone ~ lamx)\` gives
\`lambdaone = exp(lamx)\`.

See the \`binegbin.R\` file header for why Negative-Binomial components
are used instead of an observation-level random effect on \[bipois()\] –
briefly, the random-effect version fails synthetic recovery, because
with one observed pair but three latent deviates per unit the excess
deviates act as residual absorbers.

Use in a brm() call as: brm( bf(y1 \| vint(y2) ~ 1, mu ~ 1 + (1 \|
vessel), nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx), lamx ~ 1, shapes
~ 1, shapex ~ 1, nl = TRUE), family = binegbin(), stanvars =
binegbin_stanvars(), data = dat )

## Usage

``` r
binegbin()

binegbin_stanvars()

log_lik_binegbin(i, prep)

posterior_predict_binegbin(i, prep, ...)

posterior_epred_binegbin(prep)
```

## Value

A brms custom_family object.

## Details

\*\*Forced \`mu\` naming, and \`y2\` via \`vint()\`.\*\* Identical
conventions to \[bipois()\] – \`mu\` is brms's mandatory dpar name, here
bound to the shared component's rate (\`lambda_shared\`), not a mean of
either response; \`y2\` travels as supplementary integer data through
\`vint()\` because \`custom_family()\` declares a single response
column. See \[bipois()\] for the full explanation, including why the
rates are spelled \`lambdaone\`/ \`lambdatwo\` in code but written
\\\lambda_1\\/ \\\lambda_2\\ in the documentation.

\*\*Order of dpars matters for the generated Stan call.\*\* brms
generates \`target += binegbin_lpmf(Y\[n\] \| mu\[n\], lambdaone\[n\],
lambdatwo\[n\], shapes\[n\], shapex\[n\], vint1\[n\])\` – dpars in the
order declared here, then vint args. \`binegbin_stan_funs\` declares
\`binegbin_lpmf\` with exactly this signature; reordering one without
the other silently swaps which rate or dispersion governs which
component.
