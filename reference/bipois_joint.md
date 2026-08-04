# Censoring-aware joint bivariate-Poisson family for brms

Censoring-aware extension of \[bipois()\]. Models the same
trivariate-reduction bivariate-Poisson pair \`(y1, y2)\` – \`y1 =
N_shared + N1\`, \`y2 = N_shared + N2\`, with \`N_shared ~
Poisson(mu)\`, \`N1 ~ Poisson(lambdaone)\`, \`N2 ~ Poisson(lambdatwo)\`
mutually independent given their rates – but allows the first margin
(\`y1\`) to be UNOBSERVED on some rows. Each row carries two
supplementary integers via \`vint()\`: \`y2\` (the always-observed
second count) and \`y1_obs\` (a 0/1 flag marking whether \`y1\` was
observed for that row).

Stands to \[bipois()\] exactly as \[binegbin_joint()\] stands to
\[binegbin()\], and is the equidispersed special case of
\[binegbin_joint()\].

On \`y1_obs == 1\` (matched) rows the likelihood is the full joint
\[bipois()\] lpmf on \`(y1, y2)\`. On \`y1_obs == 0\` (\`y2\`-only) rows
it is the \`y1\`-integrated marginal of that same joint, which for
Poisson components is available in closed form: a sum of independent
Poissons is Poisson, so \`y2 ~ Poisson(mu + lambdatwo)\` exactly.
\[binegbin_joint()\] must evaluate the corresponding convolution as a
sum; this family does not.

Three dpars, the same as \[bipois()\]: \`mu\` (shared rate),
\`lambdaone\` and \`lambdatwo\` (the two source-specific rates), all
\`link = "log"\`.

Use in a brm() call as: brm( bf(y1 \| vint(y2, y1_obs) ~ 1, mu ~ 1 + (1
\| vessel) + (1 \| vessel:trip_id), nlf(lambdaone ~ lamx + methd),
nlf(lambdatwo ~ lamx - methd), lamx ~ 1, methd ~ 1, nl = TRUE), family =
bipois_joint(), stanvars = bipois_joint_stanvars(), data = dat )

## Usage

``` r
bipois_joint()

bipois_joint_stanvars()

log_lik_bipois_joint(i, prep)

posterior_predict_bipois_joint(i, prep, ...)

posterior_epred_bipois_joint(prep)
```

## Value

A brms custom_family object.

## Details

\*\*When to use this rather than \[binegbin_joint()\].\*\* This family
fixes each latent component's variance equal to its mean. Where the
counts are genuinely overdispersed relative to that,
\[binegbin_joint()\] is the correct model and this one will understate
the marginal variances. Where they are not, \[binegbin_joint()\] can
only represent the fit by driving its dispersions to their Poisson limit
(\`shape\` \\\to\infty\\, equivalently \`kappa\` \\\to 0\\), a boundary
at which sampling degrades; fitting the equidispersed family directly
avoids it. Compare the two with \`loo()\`.

\*\*What each rate is identified from.\*\* \`lambdatwo\` appears on both
branches, so every row informs it. \`lambdaone\` appears only on the
matched branch and is identified solely by the matched rows. \`mu\`
appears on both, but the censored branch sees it only through the sum
\`mu + lambdatwo\` – those rows constrain the total rate of the observed
margin, not how it divides between the shared and source-2-only
components. Separating \`mu\` from \`lambdatwo\` – and so estimating the
congruence \\f\\ – therefore also rests on the matched rows. With few of
them, \`mu\` and \`lambdatwo\` trade off along their sum and the prior
does correspondingly more of the work, however many censored rows the
design contains.

\*\*Two \`vint()\` arguments, in declared order.\*\* brms appends
\`vint()\` integers to the generated lpmf call in the order they are
listed in the formula's \`vint()\` term, matching the \`vars\` declared
here (\`c("vint1\[n\]", "vint2\[n\]")\`): so \`vint(y2, y1_obs)\` binds
\`vint1 = y2\` and \`vint2 = y1_obs\`. brms generates \`target +=
bipois_joint_lpmf(Y\[n\] \| mu\[n\], lambdaone\[n\], lambdatwo\[n\],
vint1\[n\], vint2\[n\])\` – dpars in the order declared here, then the
two vint args. \`bipois_joint_stan_funs\` declares \`bipois_joint_lpmf\`
with exactly this signature; reordering the dpars or the two \`vint()\`
terms without matching the Stan signature silently swaps which rate
governs which component or which integer is the branch flag.

On \`y1_obs == 0\` rows the response column \`Y\` is not read by the
likelihood, so any integer placeholder there is inert. Supply one rather
than \`NA\`, which brms rejects before the family is reached.

\*\*Forced \`mu\` naming, and the second count via \`vint()\`.\*\*
Identical conventions to \[bipois()\] – \`mu\` is brms's mandatory dpar
name, here bound to the shared component's rate, not a mean of either
response; \`y2\` (and \`y1_obs\`) travel as supplementary integer data
through \`vint()\` because \`custom_family()\` declares a single
response column. See \[bipois()\] for the full explanation, including
why the rates are spelled \`lambdaone\`/\`lambdatwo\` in code but
written \\\lambda_1\\/\\\lambda_2\\ in the documentation.
