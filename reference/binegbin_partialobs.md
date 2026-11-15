# Joint bivariate family with negative-binomial components for partially observed pairs

Returns a version of
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
for partially paired designs, in which source 1 has missing values. The
generative model, family name, dpars, and post-processing methods are
the same as for
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md).
Each row also supplies an observation indicator as the second argument
of `vint()`, after `y2`:

    brm(
      bf(y1 | vint(y2, y1_obs) ~ ...),
      family   = binegbin_partialobs(),
      stanvars = binegbin_partialobs_stanvars(),
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
binegbin_partialobs()

binegbin_partialobs_stanvars()
```

## Value

`binegbin_partialobs()` returns a brms `custom_family` object.
`binegbin_partialobs_stanvars()` returns a `stanvars` object holding the
Stan code for the corresponding `_lpmf`. The returned family has the
same `name` as
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
so the post-processing methods documented there apply to a fit made with
either constructor.

## Details

On a paired row, the likelihood term is the joint probability of
`(y1, y2)` under
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md).
On an unpaired row, the likelihood term is the marginal probability of
`y2`, the joint probability summed over every possible value of the
unrecorded `y1`:
`P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo)`
over `k = 0, ..., y2`.

To impute `y1` on the unpaired rows, use
[`posterior_predict_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md)
or
[`posterior_epred_y1_given_y2()`](https://anhsmith.github.io/bicountbrms/reference/posterior_predict_y1_given_y2.md);
see
[`vignette("prediction")`](https://anhsmith.github.io/bicountbrms/articles/prediction.md).

`lambdaone` and `shapexone` do not appear in the likelihood term for an
unpaired row, so `lambdaone` and `shapexone` are estimated from the
paired rows alone. With few paired rows, the posteriors of `lambdaone`
and `shapexone` are wide and depend strongly on their priors. `mu`,
`shapes`, `lambdatwo`, and `shapextwo` appear in both likelihood terms
and are informed by every row.

## See also

[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
for the fully paired case;
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
for the equidispersed Poisson counterpart.
