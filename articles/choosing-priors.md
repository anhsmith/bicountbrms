# Choosing priors

``` r

library(brms)
library(bicountbrms)
```

brms gives a custom family’s `mu` a default `student_t` prior and leaves
every other distributional parameter flat and improper. For a six-dpar
[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
model,
[`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
returns

                      prior     class      dpar
     student_t(3, 1.6, 2.5) Intercept
                            Intercept lambdaone
                            Intercept lambdatwo
                            Intercept    shapes
                            Intercept shapexone
                            Intercept shapextwo

An empty `prior` column means an improper uniform on the whole real
line. That is not brms behaving unusually: for its built-in families the
remaining parameters are things like a residual standard deviation, for
which brms does supply a default. A custom family has no such defaults,
because brms cannot know what its parameters mean. A reader who assumes
otherwise, reasonably enough from experience with
[`gaussian()`](https://rdrr.io/r/stats/family.html) or
[`negbinomial()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
has no reason to run
[`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
at all. That assumption is why this page exists.

## Every parameter identified by the matched rows alone is left improper

Under partial observation the two constructors do not identify every
parameter equally well. `lambdaone` and `shapexone` appear only on the
matched branch of the likelihood, as does any between-source bias, so
the rows on which the first count was recorded are the only rows that
inform them. `mu`, `shapes`, `lambdatwo` and `shapextwo` appear on both
branches, and every row informs them.

Set that division beside the
[`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
output above. Every parameter in the first group is left without a
prior. `shapes` and `shapextwo` are left without one as well, although
every row informs them, so the containment runs one way only. The one
parameter brms does supply a prior for, `mu`, is informed by every row.
Under a design in which few rows are matched, the parameters whose
posteriors depend most on their priors are all parameters for which no
prior has been set.

For
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
the statement needs one qualification. On an unmatched row that family’s
likelihood involves `mu` only through the sum `mu + lambdatwo`, so those
rows constrain the total rate of the observed margin without dividing it
between the shared and the source-specific component. Separating `mu`
from `lambdatwo`, and so estimating the congruence , is also informed by
the matched rows alone.

## The symptom is divergent transitions

Left unregularised, the chains explore the flat tail of an improper
prior on a log-linked positive parameter, where the likelihood is nearly
level and the posterior therefore follows the prior. The
asymmetric-dispersion fit in
`tests/testthat/test-binegbin-dispersions.R` required priors on the two
rates and the three dispersions: run without them it produced 680
divergent transitions and a maximum of 1.54, which is to say chains that
never mixed.

Divergences are the diagnostic to watch here rather than alone. This
package applies an gate of 1.02 rather than the 1.01 of Vehtari et al.
([2021](#ref-vehtari2021)), for reasons set out in
`tests/testthat/helper-coverage.R`, and a model can sit inside that gate
while still producing divergent transitions.

## A weakly informative prior on each log-scale parameter

Every distributional parameter of every constructor in this package is
log-linked, so a prior written on the log scale is lognormal on the
natural scale. Check its quantiles on the natural scale, where the
parameter is interpreted ([McElreath 2020](#ref-mcelreath2020)).

``` r

prior = c(
  prior(normal(2, 1), class = "Intercept"),                     # log shared rate
  prior(normal(2, 1), class = "Intercept", dpar = "lambdaone"),
  prior(normal(2, 1), class = "Intercept", dpar = "lambdatwo"),
  prior(normal(2, 1), class = "Intercept", dpar = "shapes"),
  prior(normal(2, 1), class = "Intercept", dpar = "shapexone"),
  prior(normal(2, 1), class = "Intercept", dpar = "shapextwo")
)
```

On the natural scale `normal(2, 1)` places 95% of its mass in roughly ,
which suits rates and dispersions of order 1 to 50.
[`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
and
[`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md)
have no dispersion parameters, so only the three rate lines apply to
them.

To adapt these priors to counts on a different scale, shift the **mean**
rather than widening the standard deviation. Because the link is
logarithmic, widening a normal prior produces a heavier-tailed
lognormal, which places *more* mass in the flat region the sampler
explores rather than less ([Smith et al.
2020](#ref-smithInstantaneousVsNoninstantaneous2020), supplement 3). A
prior centred at the scale of the data and no wider than it needs to be
is both more informative about what is plausible and better behaved.

## The class, dpar and nlpar slots

The slot a prior occupies depends on how the parameter reached the
model, not on what the parameter means. A dispersion supplied directly
as a distributional parameter takes `class = "Intercept"` and `dpar`. A
rate supplied through
[`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
takes `class = "b"` and `nlpar`.

``` r

# lambdaone and lambdatwo tied to one value through a non-linear parameter
bf(y1 | vint(y2) ~ 1, mu ~ 1,
   nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx),
   lamx ~ 1, shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE)

prior = c(
  prior(normal(2, 1), class = "Intercept"),                     # mu
  prior(normal(2, 1), class = "b",         nlpar = "lamx"),     # not dpar
  prior(normal(2, 1), class = "Intercept", dpar  = "shapes")
)
```

[`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
on the formula returns the combination each parameter requires, and is
the reliable way to find it. A prior naming a slot the model does not
contain is dropped without a warning, which leaves the parameter
improper and the symptom above unexplained.

That failure mode has one especially easy instance. To impose the two
excess dispersions are routed through a single non-linear parameter,
described in [The families and their
parameters](https://anhsmith.github.io/bicountbrms/articles/families-and-parameters.md):

``` r

bf(y1 | vint(y2) ~ 1, mu ~ 1,
   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
   lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE)

prior(normal(0, 1.5), class = "b", nlpar = "shapexx")
```

Both fields differ from the pre-0.10.0 spelling, which was
`class = "Intercept", dpar = "shapex"`. A prior written the old way
names no parameter in the tied model, is dropped silently, and leaves
`shapexx` flat.

## Shrinkage priors in the (M, f, δ) coordinates

The three rates can be reparameterised as overall level , congruence and
source bias , with each dispersion on the standard deviation scale . In
those coordinates every parameter has a finite null: is no bias between
the two sources, and is the Poisson limit exactly, rather than the of
the native scale.

A shrinkage prior requires a finite null. A half-normal on and a normal
on place their maximum at the simpler model and decay away from it, so
the data rather than the prior must move the fit. They follow the
principles of Simpson et al.
([2017](#ref-simpsonPenalisingModelComponent2017)), a base model at zero
and a density decaying away from it, without being penalised-complexity
priors in the strict sense.

[The anatomy of a paired
count](https://anhsmith.github.io/bicountbrms/articles/paired-count-anatomy.md)
sets out that argument in full, shows the prior pushforward onto the
native scale, and fits a model in those coordinates. Read it after this
page rather than instead of it: the recipes above apply whichever
parameterisation is used, because a prior in still occupies one of the
three slots described here.

## References

McElreath, Richard. 2020. *Statistical Rethinking: A Bayesian Course
with Examples in R and Stan*. 2nd ed. Chapman; Hall/CRC.
<https://doi.org/10.1201/9780429029608>.

Simpson, Daniel, Håvard Rue, Andrea Riebler, Thiago G. Martins, and
Sigrunn H. Sørbye. 2017. “Penalising Model Component Complexity: A
Principled, Practical Approach to Constructing Priors.” *Statistical
Science* 32 (1): 1–28. <https://doi.org/10.1214/16-STS576>.

Smith, A. N. H., D. Acuña-Marrero, P. Salinas-de-León, E. S. Harvey, M.
D. M. Pawley, and M. J. Anderson. 2020. “Instantaneous Vs.
Non-Instantaneous Diver-Operated Stereo-Video (DOV) Surveys of Highly
Mobile Sharks in the Galápagos Marine Reserve.” *Marine Ecology Progress
Series* 649: 111–23. <https://doi.org/10.3354/meps13447>.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $`\widehat{R}`$ for Assessing Convergence of
MCMC (with Discussion).” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.
