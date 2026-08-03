# Censoring-aware joint bivariate-Negative-Binomial family for brms

Censoring-aware extension of \[binegbin()\]. Models the same trivariate-
reduction bivariate Negative-Binomial pair \`(y1, y2)\` – \`y1 =
N_shared + N1\`, \`y2 = N_shared + N2\`, with \`N_shared ~ NB2(mu,
shapes)\`, \`N1 ~ NB2(lambdaone, shapexone)\`, \`N2 ~ NB2(lambdatwo,
shapextwo)\` mutually independent given their rates – but allows the
first margin (\`y1\`) to be UNOBSERVED on some rows. Each row carries
two supplementary integers via \`vint()\`: \`y2\` (the always-observed
second count) and \`y1_obs\` (a 0/1 flag marking whether \`y1\` was
observed for that row).

On \`y1_obs == 1\` (matched) rows the likelihood is the full joint
\[binegbin()\] lpmf on \`(y1, y2)\`. On \`y1_obs == 0\` (\`y2\`-only)
rows it is the \`y1\`-integrated marginal of that same joint, \`P(y2) =
sum_k NB2(k \| mu, shapes) NB2(y2 - k \| lambdatwo, shapextwo)\` – NOT a
separate single-dispersion \`neg_binomial_2\` on \`y2\`, which would be
a different model inconsistent with the matched decomposition. This lets
one \`brm()\` call pool matched and \`y2\`-only rows under one coherent
likelihood: \`lambdaone\`, the between-source bias and \`shapexone\` are
identified only by the matched rows, while the \`y2\`-only rows sharpen
\`mu\`, \`shapes\`, \`lambdatwo\`, \`shapextwo\`, and the shared
vessel/trip random-effect structure.

Six dpars: the three rates (\`mu\` = shared rate,
\`lambdaone\`/\`lambdatwo\` = the two source-specific rates) plus three
dispersions – \`shapes\` for the shared component and
\`shapexone\`/\`shapextwo\` for the two source-specific excess
components. All six use \`link = "log"\` (see \[binegbin()\]). To share
a level across a pair of dpars and split them by a directional bias,
supply them through non-linear formulas \*without\* an explicit
\`exp()\` – the log link applies it, so \`nlf(lambdaone ~ lamx +
methd)\` gives \`lambdaone = exp(lamx + methd)\`.

Use in a brm() call as: brm( bf(y1 \| vint(y2, y1_obs) ~ 1, mu ~ 1 + (1
\| vessel) + (1 \| vessel:trip_id), nlf(lambdaone ~ lamx + methd),
nlf(lambdatwo ~ lamx - methd), lamx ~ 1, methd ~ 1, shapes ~ 1,
shapexone ~ 1, shapextwo ~ 1, nl = TRUE), family = binegbin_joint(),
stanvars = binegbin_joint_stanvars(), data = dat )

## Usage

``` r
binegbin_joint()

binegbin_joint_stanvars()

log_lik_binegbin_joint(i, prep)

posterior_predict_binegbin_joint(i, prep, ...)
```

## Value

A brms custom_family object.

## Details

\*\*The symmetric model is a formula constraint.\*\* Before 0.8.0 this
family carried a single excess dispersion \`shapex\` shared by both
margins. That model is the constraint \`shapexone == shapextwo\`,
obtained by routing both through one non-linear parameter:

“\` bf(y1 \| vint(y2, y1_obs) ~ 1, mu ~ 1 + (1 \| vessel), nlf(shapexone
~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1, ..., nl = TRUE) “\`

The resulting likelihood is term-for-term the pre-0.8.0 five-dpar one (a
package test pins this). Stored five-dpar fits keep working unchanged:
their single \`shapex\` resolves to both \`shapexone\` and \`shapextwo\`
when post-processed. See \`migration/family-unification.md\` for the
migration detail.

\*\*What each dispersion is identified from.\*\* \`shapextwo\` governs
the always-observed margin and so appears on both branches; the
\`y2\`-only rows inform it. \`shapexone\` appears only on the matched
branch and is identified solely by the matched rows, as \`lambdaone\`
is. With few matched rows the posterior for \`shapexone\` is
correspondingly dominated by its prior.

\*\*Two \`vint()\` arguments, in declared order.\*\* brms appends
\`vint()\` integers to the generated lpmf call in the order they are
listed in the formula's \`vint()\` term, matching the \`vars\` declared
here (\`c("vint1\[n\]", "vint2\[n\]")\`): so \`vint(y2, y1_obs)\` binds
\`vint1 = y2\` and \`vint2 = y1_obs\`. brms generates \`target +=
binegbin_joint_lpmf(Y\[n\] \| mu\[n\], lambdaone\[n\], lambdatwo\[n\],
shapes\[n\], shapexone\[n\], shapextwo\[n\], vint1\[n\], vint2\[n\])\` –
dpars in the order declared here, then the two vint args.
\`binegbin_joint_stan_funs\` declares \`binegbin_joint_lpmf\` with
exactly this signature; reordering the dpars or the two \`vint()\` terms
without matching the Stan signature silently swaps which quantity
governs which component or which integer is the branch flag.

\*\*Forced \`mu\` naming, and the second count via \`vint()\`.\*\*
Identical conventions to \[binegbin()\]/\[bipois()\] – \`mu\` is brms's
mandatory dpar name, here bound to the shared component's rate, not a
mean of either response; \`y2\` (and \`y1_obs\`) travel as supplementary
integer data through \`vint()\` because \`custom_family()\` declares a
single response column. See \[bipois()\] for the full explanation,
including why the rates are spelled \`lambdaone\`/\`lambdatwo\` in code
but written \\\lambda_1\\/\\\lambda_2\\ in the documentation.

\*\*Relationship to \[binegbin()\].\*\* On \`y1_obs == 1\` rows with
\`shapexone == shapextwo\` this family's lpmf equals the \[binegbin()\]
lpmf exactly (same marginalisation sum). The \`y1_obs == 0\` branch is
the \`y1\`-integrated marginal of that same bivariate model. The package
tests pin both identities (marginal identity; binegbin equivalence).
