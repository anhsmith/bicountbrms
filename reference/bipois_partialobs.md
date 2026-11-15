# Joint bivariate-Poisson family for partially observed pairs

Returns a version of
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
for partially paired designs, in which source 1 has missing values. The
generative model, family name, dpars, and post-processing methods are
the same as for
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md).
Each row also supplies an observation indicator as the second argument
of `vint()`, after `y2`:

    brm(
      bf(y1 | vint(y2, y1_obs) ~ ...),
      family   = bipois_partialobs(),
      stanvars = bipois_partialobs_stanvars(),
      data     = dat
    )

This template is for the native parameterisation. To fit the model in
the congruence–bias parameterisation, which we recommend, write the
formula with
[`bf_mfd()`](https://anhsmith.github.io/bicountbrms/reference/bf_mfd.md).

`y1_obs` is `1` on rows where both counts were recorded and `0` on rows
where only `y2` was recorded. An indicator other than `0` or `1` makes
the sampler fail at initialisation, with a message in the chain output
that names the cause. The check also catches the two `vint()` arguments
written in reverse order.

**Do not use `NA` for missing `y1` values.** brms drops any row with a
missing response with only a warning. On a row with `y1_obs == 0`, the
likelihood term does not depend on `y1`, so `y1` can be any integer. We
suggest `-1`, because it cannot be confused with a recorded count:

    dat$y1[dat$y1_obs == 0] <- -1L

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

On a paired row, the likelihood term is the joint probability of
`(y1, y2)` under
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md).
On an unpaired row, the likelihood term is the marginal probability of
`y2`, the joint probability summed over every possible value of the
unrecorded `y1`. That sum has the closed form
`y2 ~ Poisson(mu + lambdatwo)`, because a sum of independent Poisson
variables is Poisson.

To impute `y1` on the unpaired rows, use
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
or
[`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md);
see
[`vignette("prediction")`](https://anhsmith.github.io/bicountbrms/articles/prediction.md).

An unpaired row informs `mu` and `lambdatwo` only through their sum,
`mu + lambdatwo`, which is the expected value of `y2`. Two quantities
are therefore estimated from the paired rows alone: (i) `lambdaone`,
which does not appear in the likelihood term for an unpaired row, and
(ii) the division of `mu + lambdatwo` between the shared and
source-2-exclusive components. The congruence \\f\\ depends on both, so
\\f\\ is also estimated from the paired rows alone. With few paired
rows, the posteriors of these quantities stay close to their priors,
however many unpaired rows the design contains.

`bipois_partialobs()` fixes the variance of each latent component equal
to the mean of that component. Where the counts are overdispersed
relative to a Poisson distribution,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
is the correct model and `bipois_partialobs()` will understate the
marginal variances. Where the counts are not overdispersed,
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
fits the counts only with `shapes`, `shapexone`, and `shapextwo` very
large (each `kappa` near 0), and sampling is poor when the dispersions
are that large. `bipois_partialobs()` has no dispersions, so it avoids
the problem. Compare the fits from the two constructors with
[`loo()`](https://mc-stan.org/loo/reference/loo.html), which scores `y1`
and `y2` on a paired row and `y2` on an unpaired row.

## See also

[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
for the fully paired case;
[`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md)
for the overdispersed counterpart.
