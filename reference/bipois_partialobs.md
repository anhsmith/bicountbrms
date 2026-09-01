# Joint bivariate-Poisson family for partially observed pairs

[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
for a design in which the first count is missing on some rows. Same
generative model, same `name`, same three dpars, same likelihood and the
same post-processing methods – the only difference is that each row
supplies a second integer, an observation flag, through `vint()`:

    bf(y1 | vint(y2, y1_obs) ~ ...)

`y1_obs` is a 0/1 integer column: `1` where both counts were recorded,
`0` where the first was not. `y1` may hold any non-negative integer on
those rows – `0` is the conventional placeholder – because the
likelihood does not read it. Do not use `NA`, which brms drops before
fitting, taking the row's observed `y2` with it.

**Contribution of a matched row and of an unmatched row.** A matched row
(`y1_obs == 1`) uses the full joint
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
lpmf on `(y1, y2)`. A row whose first count was never recorded
(`y1_obs == 0`) contributes the second count's marginal *from the same
model*. For Poisson components that marginal is closed form – a sum of
independent Poissons is Poisson – so it is exactly
`y2 ~ Poisson(mu + lambdatwo)`.
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
must evaluate the corresponding convolution as a sum; this family does
not. Either way the row is not dropped and is not given a different
model: it still informs `mu`, `lambdatwo` and any group-level effects.

**Imputation after fitting.** The fitted model can impute the unobserved
first count conditional on the observed second one, which is usually why
someone wanted this.
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
return a `y1` draw and `E[y1 | y2]` for *every* row, matched and
unmatched alike – `y1_obs` selects a likelihood branch, not a
prediction.

**A design consequence, worth knowing before the data are collected.**
`lambdatwo` appears on both branches, so every row informs it.
`lambdaone` appears only on the matched branch and is identified by the
matched rows *alone*. `mu` appears on both, but the unmatched branch
sees it only through the sum `mu + lambdatwo` – those rows constrain the
total rate of the observed margin, not how it divides between the shared
and source-2-only components. Separating `mu` from `lambdatwo`, and so
estimating the congruence \\f\\, is therefore also informed by the
matched rows. With few of them, `mu` and `lambdatwo` trade off along
their sum and the prior does correspondingly more of the work, however
many unmatched rows the design contains.

**This is not censoring in brms's sense.** brms's `cens()` addition term
means a value known to lie in a set – `left`, `right`, `interval`. Here
the first count is not observed at all and the likelihood marginalises
over its whole support. This family was called `bipois_cens()` up to
0.9.1; the name was wrong and was changed at 0.10.0. Do not combine this
family with `cens()`.

Use in a brm() call as: brm( bf(y1 \| vint(y2, y1_obs) ~ 1, mu ~ 1 + (1
\| vessel) + (1 \| vessel:trip_id), nlf(lambdaone ~ lamx + methd),
nlf(lambdatwo ~ lamx - methd), lamx ~ 1, methd ~ 1, nl = TRUE), family =
bipois_partialobs(), stanvars = bipois_partialobs_stanvars(), data = dat
)

## Usage

``` r
bipois_partialobs()

bipois_partialobs_stanvars()
```

## Value

`bipois_partialobs()` returns a brms `custom_family` object.
`bipois_partialobs_stanvars()` returns a `stanvars` object holding the
Stan code for the corresponding `_lpmf`. The returned family has the
same `name` as
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
so the post-processing methods documented there apply to a fit made with
either constructor.

## Details

**Choosing between this family and
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).**
This family fixes each latent component's variance equal to its mean.
Where the counts are genuinely overdispersed relative to that,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
is the correct model and this one will understate the marginal
variances. Where they are not,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
can only represent the fit by driving its dispersions to their Poisson
limit (`shape` \\\to\infty\\, equivalently `kappa` \\\to 0\\), a
boundary at which sampling degrades; fitting the equidispersed family
directly avoids it. Compare the two with
[`loo()`](https://mc-stan.org/loo/reference/loo.html).

**One likelihood, two constructors.** This returns the same
`custom_family` `name` as
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
so brms resolves both to one `bipois_lpmf` and one set of
[`log_lik_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
/
[`posterior_predict_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
/
[`posterior_epred_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
methods. The matched branch is therefore not a second copy that can
drift from the fully paired likelihood; it is the same code. See
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
for why `vars` declares a literal in the plain constructor rather than
the two families declaring overloaded Stan functions.

**Two `vint()` arguments, in declared order.** brms appends `vint()`
integers to the generated lpmf call in the order they are listed in the
formula's `vint()` term, matching the `vars` declared here
(`c("vint1[n]", "vint2[n]")`): so `vint(y2, y1_obs)` binds `vint1 = y2`
and `vint2 = y1_obs`. brms generates
`target += bipois_lpmf(Y[n] | mu[n], lambdaone[n], lambdatwo[n], vint1[n], vint2[n])`.
Reordering the dpars or the two `vint()` terms without matching the Stan
signature silently swaps which rate governs which component or which
integer is the branch flag.

**Distinguishing the two constructors in a stored fit.** `family$name`
is `"bipois"` either way. What distinguishes them is the presence of the
second supplementary integer:

    "vint2" %in% names(brms::standata(fit))   # TRUE for a partially observed fit
    fit$family$vars     # c("vint1[n]", "vint2[n]") or c("vint1[n]", "1")

## See also

[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
for the fully paired case;
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
for the overdispersed counterpart.
