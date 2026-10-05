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
fitting, taking the observed `y2` on that row with it.

**Contribution of a paired row and of an unpaired row.** A paired row
(`y1_obs == 1`) uses the full joint
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
lpmf on `(y1, y2)`. A row whose first count was never recorded
(`y1_obs == 0`) contributes the marginal of the second count *from the
same model*. For Poisson components that marginal is closed form – a sum
of independent Poissons is Poisson – so it is exactly
`y2 ~ Poisson(mu + lambdatwo)`.
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
evaluates the marginal of `y2`, the convolution of the shared and
source-2-exclusive negative-binomial components, as a finite sum over
`k = 0..y2`. In both families, an unpaired row is not dropped and is not
given a different model: it still informs `mu`, `lambdatwo` and any
group-level effects.

**Imputation after fitting.** The fitted model can impute the unobserved
first count conditional on the observed second count.
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
return a `y1` draw and `E[y1 | y2]` for *every* row, paired and unpaired
alike – `y1_obs` selects a likelihood branch, not a prediction.

**Parameters identified only by the paired rows.** `lambdatwo` appears
on both branches, so every row informs `lambdatwo`. `lambdaone` appears
only on the paired branch and is identified by the paired rows alone.
The likelihood term for an unpaired row depends on `mu` and `lambdatwo`
only through `mu + lambdatwo`. Unpaired rows therefore constrain the
expected count of the second source, but not how `mu + lambdatwo`
divides between the shared and source-2-exclusive components. The
separation of `mu` from `lambdatwo`, and so the estimate of the
congruence \\f\\, is informed by the paired rows alone. With few paired
rows, the likelihood is nearly flat along `mu + lambdatwo` = constant,
so in that direction the posterior follows the prior, however many
unpaired rows the design contains.

**Relation to censoring in brms.** brms uses the `cens()` addition term
for a value known to lie in a set (`left`, `right` or `interval`). On an
unpaired row, the first count is not observed at all, and the likelihood
marginalises over its whole support. `bipois_partialobs()` was called
`bipois_cens()` up to 0.9.1. Do not combine `bipois_partialobs()` with
`cens()`.

Use in a brm() call as: brm( bf(y1 \| vint(y2, y1_obs) ~ 1, mu ~ 1 + (1
\| vessel) + (1 \| vessel:trip_id), nlf(lambdaone ~ lamx + delta),
nlf(lambdatwo ~ lamx - delta), lamx ~ 1, delta ~ 1, nl = TRUE), family =
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

**Choosing between `bipois_partialobs()` and
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md).**
`bipois_partialobs()` fixes the variance of each latent component equal
to its mean. Where the counts are overdispersed relative to a Poisson
distribution,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
is the correct model and `bipois_partialobs()` will understate the
marginal variances. Where the counts are not overdispersed,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
fits such counts only by driving its dispersions to their Poisson limit
(`shape` \\\to\infty\\, equivalently `kappa` \\\to 0\\), a boundary at
which sampling degrades; fitting `bipois_partialobs()` directly avoids
that boundary. Compare the fits from the two constructors with
[`loo()`](https://mc-stan.org/loo/reference/loo.html).

**Shared Stan function and post-processing methods.**
`bipois_partialobs()` returns the same `custom_family` `name` as
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md),
so brms resolves both constructors to one `bipois_lpmf` and one set of
[`log_lik_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
/
[`posterior_predict_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
/
[`posterior_epred_bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
methods. The paired branch of `bipois_partialobs()` therefore runs the
same code as
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md).
See
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
for why `vars` declares a literal in the plain constructor rather than
the two constructors declaring overloaded Stan functions.

**Two `vint()` arguments, in declared order.** brms appends `vint()`
integers to the generated lpmf call in the order they are listed in the
`vint()` term of the formula, matching the `vars` declared here
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
