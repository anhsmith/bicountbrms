# Changelog

## pairedcountbrms 0.8.0

- **Breaking (dpar split):
  [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)’s
  single excess dispersion `shapex` is now the per-margin pair
  `shapexone`/`shapextwo`**, giving the family six dpars (`mu`,
  `lambdaone`, `lambdatwo`, `shapes`, `shapexone`, `shapextwo`), all
  log-linked. The two source-only excess components may now differ in
  overdispersion, which the single `shapex` forced them not to. The
  rates were already free to differ; the dispersions now are too.

- **The symmetric model is a formula constraint, not a separate
  family.** `shapexone == shapextwo` recovers the pre-0.8.0 likelihood
  term for term. Express it by routing both dpars through one non-linear
  parameter:

  ``` r

  bf(y1 | vint(y2, y1_obs) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(shapexone ~ shapexx),
     nlf(shapextwo ~ shapexx),
     shapexx ~ 1, ..., nl = TRUE)
  ```

  This is the only source change a symmetric model needs. Note the prior
  moves with the name: `prior(..., dpar = "shapex")` becomes
  `prior(..., nlpar = "shapexx")`.

- **Fits made before this release still post-process; no refitting and
  no shim are required.** brms resolves a custom family’s `log_lik_*` /
  `posterior_predict_*` by name against the live search path at call
  time, not from anything frozen in the fit, so a stored fit always runs
  the currently attached code. Only the dpar NAMES it carries are
  frozen, and both excess dispersions now fall back to a five-dpar fit’s
  single `shapex` — which is precisely the constraint that fit was
  estimated under. Verified on stored five-dpar fits:
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  and
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  agree with 0.7.0 bitwise, and
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) runs unchanged.

- The generalised family absorbs a project-local asymmetric family
  developed outside the package. During the migration the two were
  compared element by element over 15,480,000 pointwise log-likelihood
  values from ten stored fits, and every value was bitwise identical.
  That was a one-off check against artefacts outside this repository,
  not a package test; it is recorded with its method and versions in
  `migration/family-unification.md`, alongside the dpar mapping table
  and adoption steps.

- [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  gains `kappaxone`/`kappaxtwo`, and
  [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_dpars_to_mfd.md)
  gains `shapexone`/`shapextwo`, so the converters serve
  [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)’s
  per-margin excess dispersions under the family’s own dpar names
  instead of leaving the caller to re-derive `shape = 1/kappa^2`. Purely
  additive: the existing `kappax`/`shapex` arguments are unchanged and
  still return `shapex`/`kappax`. Supplying `kappax` together with
  either of the new arguments is an error, since they are two spellings
  of the same quantity for different families.

- Note on identifiability: `shapextwo` governs the always-observed
  margin and appears on both branches, so the censored rows inform it.
  `shapexone` appears only on the matched branch and is identified
  solely by the matched rows, as `lambdaone` is. A design with few
  matched rows will learn `shapexone` weakly.

## pairedcountbrms 0.7.0

- **Breaking (dpar rename): the joint families’ excess-rate dpars
  `lambdaem`/`lambdalb` are now `lambdaone`/`lambdatwo`**, affecting
  [`bipois()`](https://anhsmith.github.io/pairedcountbrms/reference/bipois.md),
  [`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  and
  [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md).
  The former names referred to electronic monitoring and vessel
  logbooks, the two sources in the project from which the package was
  extracted, and carried no meaning outside it. The new names are
  positional: source 1 is the response, source 2 is supplied via
  `vint()`.

  They are spelled out rather than written `lambda1`/`lambda2` because
  [`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
  rejects dpar names ending in a digit, as well as dots and underscores.
  The documentation writes them as $`\lambda_1`$ and $`\lambda_2`$; a
  notation table in the README gives the correspondence.

- **Fits made before this release still post-process; no refitting is
  required.** A `brmsfit` stores its own family object, so
  [`prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html)
  on an older fit returns a `prep` carrying the previous dpar names.
  Rate reads in the joint families now resolve the name against the fit,
  so
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html),
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  return numerically identical results on pre-0.7.0 fits
  (`tests/testthat/test-dpar-compat.R`).

  Code that names the dpars does require editing:
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)/[`lf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  terms and `prior(..., dpar = )` arguments. Data column names are
  unaffected, since the second count and the observation flag reach the
  family positionally through `vint()`.

- **Breaking (argument names):**
  [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_dpars_to_mfd.md)
  takes `lambdaone`/`lambdatwo`, and
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  returns a list with those names. Positional calls are unaffected.

- **Documentation notation unified.** The joint families previously used
  $`y_{\mathrm{em}}, y_{\mathrm{lb}}, N_{10}, N_{01}`$ while the
  difference families used $`y_1, y_2`$, so the identity between the
  bivariate Poisson’s difference and the Skellam distribution was stated
  in one notation and demonstrated in another. All documentation now
  uses $`y_1, y_2`$, $`N_1, N_2`$ and $`\lambda_1, \lambda_2`$, making
  $`d = y_1 - y_2 = N_1 - N_2`$ continuous across both suites.

- **Recovery tests separated into smoke gates and a calibration
  assessment.** The recovery tests asserted that the true value fell
  within a single fit’s 90% credible interval. For a correct model with
  a calibrated posterior this is a Bernoulli(0.9) draw, failing 10% of
  the time by construction; across the ~18 such assertions in the suite,
  spurious failures were the norm. Two assertions in `test-binegbin.R`
  had been failing on every run since they were written, identically
  before and after the dpar rename, undetected because `skip_on_cran()`
  excludes the fitting tests from `R CMD check`.

  Single-fit checks are now smoke gates: convergence plus a wide (99%)
  interval, which detects gross mis-specification without failing a
  correct model by chance. Calibration is assessed separately, as the
  proportion of repeated simulate-and-refit replicates whose nominal
  interval contains the truth, with the pass threshold derived from the
  Binomial(*R*, 0.9) null (0.16% false-failure rate at *R* = 10, against
  10% per assertion previously).

  The coverage assessment costs minutes and is opt-in via
  `PAIREDCOUNTBRMS_COVERAGE=true`, a named gate rather than a silent
  default. It compiles the model once and reuses it through
  `update(recompile = FALSE)`, which keeps 10 replicates to
  approximately 3.5 minutes. Shared helpers are in
  `tests/testthat/helper-coverage.R`.

- **`binegbin` recovery and group-level composition are now separate
  tests, with less extreme generative truths.** A single test previously
  did both, on data with a vessel effect over 8 levels and severe
  overdispersion (the shared component’s variance was 40 against a mean
  of 8). Three variance channels — the group-level SD, `shapes` and
  `shapex` — competed for the same residual. `shapex` is identified only
  through the difference’s variance, and with the excess rate also
  unknown the two trade off along `shapex` = `lambda`^2 / (*V* −
  `lambda`) at fixed *V* = Var(*d*)/2, so a 1.5 SD fluctuation in
  Var(*d*) displaced the posterior substantially.

  Recovery is now asserted on an intercepts-only fit at *n* = 400 with
  approximately half of each component’s variance from overdispersion
  (`shapes` = 8, `shapex` = 3), identifiable from both directions. A
  second test adds a group-level term and checks only that the family
  composes with brms’s random-effects machinery and samples cleanly,
  making no dispersion-recovery claim. Recovery also includes a
  posterior-predictive check on Var(*d*), the quantity `shapex` governs.

- **References to the originating client project removed.** The
  documentation pointed readers at that project by name and at specific
  internal `.qmd` files — 26 references across the R sources, four help
  pages and three test files, several of which rendered onto the public
  documentation site, and none accessible to a reader. Where a pointer
  carried an argument (the regression-to-the-mean rationale for
  modelling the pair jointly, the observation-level random-effect
  failure behind the Negative-Binomial components, the prior-scale
  translation on
  [`skellam1()`](https://anhsmith.github.io/pairedcountbrms/reference/skellam1.md)),
  that argument is now stated inline. Validation claims cite the
  package’s own test files.

- **Continuous integration now runs the tests.** The repository
  previously contained one workflow, `pkgdown.yaml`, which built and
  deployed the documentation site and ran no tests. `R-CMD-check.yaml`
  adds two jobs on different cadences: `check` on every push and pull
  request, without a Stan backend (~5 minutes), and `check-stan` weekly
  and on demand, which sets `NOT_CRAN=true` and installs both rstan and
  cmdstanr so that the model fits execute.

  The fitting tests are both the most likely to detect a regression and
  too slow for every push. `skip_on_cran()` correctly excludes them from
  a routine check, but with nothing scheduled to lift it they had never
  run in continuous integration at all.

- **The $`(M, f, \delta)`$ reparameterisation is now demonstrated rather
  than described.** The anatomy article previously carried an
  unevaluated
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  fragment in place of a worked example — no `family`, no `stanvars`, no
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call, and
  never executed. It now simulates from known coordinates, fits,
  recovers all five, and converts the posterior back through
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md).

  The fit is parameterised in $`\kappa = 1/\sqrt{\phi}`$ rather than
  $`\phi`$, via `nlf(shapes ~ -2 * log(kappas))`. In $`\phi`$ the
  Poisson limit is $`\phi \to \infty`$, so there is no finite point to
  shrink towards and an exponential prior on $`\phi`$ shrinks towards
  maximum overdispersion. In $`\kappa`$ the limit is zero, which admits
  penalised-complexity priors (Simpson et al. 2017):
  [`exponential()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  on each $`\kappa`$ shrinking to Poisson, and `double_exponential()` on
  $`\delta`$ shrinking to no between-source bias. Congruence $`f`$ is
  deliberately not shrunk in either direction, both of its ends being
  degenerate.

  Each prior is plotted on the scale it is interpreted on rather than
  the linear-predictor scale it is stated on.

- **The anatomy article is now precompiled**, as the vignette is,
  because it fits a model.
  `vignettes/articles/paired-count-anatomy.Rmd.orig` is the source; the
  `.Rmd` and `vignettes/articles/figure/` are generated. Edit the
  `.orig`. `vignettes/precompile.R` handles both documents, preserving
  the property that pkgdown requires no Stan backend.

- **New README section on priors.**
  [`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
  on a
  [`binegbin()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md)
  model shows that brms assigns the `mu` dpar a default `student_t`
  prior but leaves `lamx`, `shapes` and `shapex` flat and improper;
  these are the weakly-identified parameters. In this package’s recovery
  test, omitting priors produced a divergent transition and an Rhat of
  1.0101. The section gives a weakly-informative set, notes that the
  `class`/`dpar`/`nlpar` slots differ between a rate supplied through
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  and the dispersions, and advises shifting the prior mean rather than
  increasing its SD. The recovery tests use these priors.

- **New README note on `nl = TRUE`.** In a non-linear brms formula the
  main formula’s right-hand side is a non-linear expression for `mu`,
  not a request for an intercept.
  `bf(y1 | vint(y2) ~ 1, nlf(...), ..., nl = TRUE)` without a `mu ~ 1`
  term therefore generates `mu[n] = exp(1)`, fixing the shared rate at
  *e*. This does not error; sampling becomes very slow. All the joint
  families’ documented examples use `nl = TRUE` so that
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  can tie the two excess rates.

- **`citation("pairedcountbrms")` now returns a usable citation.** There
  was no `inst/CITATION`, and the auto-generated fallback could not
  determine a year from `DESCRIPTION`, rendering as `Smith A (????)`
  with a warning. `inst/CITATION` reads the version and year from
  `DESCRIPTION`, and a `Date` field has been added. The header directs
  users to `citation("brms")` and to the distribution papers in the
  README.

- Examples and vignettes use `y1`/`y2` for the two count columns and
  `y1_obs` for
  [`binegbin_joint()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_joint.md)’s
  observation flag, previously `y_em`/`y_lb`/`em_obs`. These are
  illustrative data column names, not API.

## pairedcountbrms 0.6.0

- **Renamed: `skellambrms` is now `pairedcountbrms`.** The old name
  named one family; the package’s subject is the comparison of two
  paired count sources, by two complementary routes — the difference
  families (`skellam1`/`skellam2`, `dnorm1`/`dnorm2`,
  `dlaplace1`/`dlaplace2`) and the joint bivariate families (`bipois`,
  `binegbin`, `binegbin_joint`). Skellam is one member of that set, and
  no longer the most used one.

  **No family names change**, so no fitted model needs refitting:
  `binegbin`, `binegbin_joint`, `bipois`, `skellam1`/`skellam2`,
  `dnorm1`/`dnorm2` and `dlaplace1`/`dlaplace2` all keep their names,
  and brms continues to resolve each fit’s `log_lik_*` /
  `posterior_predict_*` / `posterior_epred_*` methods off the attached
  search path exactly as before. The only change a user needs to make is
  [`library(skellambrms)`](https://rdrr.io/r/base/library.html) →
  [`library(pairedcountbrms)`](https://github.com/anhsmith/pairedcountbrms)
  (and any `skellambrms::` prefix).

  The GitHub repository moves to `anhsmith/pairedcountbrms`. GitHub
  serves a permanent redirect from the old path for both web and git, so
  existing clones and `pak::pak("anhsmith/skellambrms")` calls keep
  working.

- New coordinate helpers,
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)
  and
  [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_dpars_to_mfd.md),
  converting between the rate dpars the joint families take (`mu`,
  `lambdaem`, `lambdalb`, plus `shapes`/`shapex`) and the interpretable
  `(M, f, delta)` coordinates — overall level, congruence, and method
  bias — along with the SD-scale dispersions `kappas`/`kappax`. Pure
  transforms; they fit nothing, and fitting in these coordinates still
  goes through
  [`nlf()`](https://paulbuerkner.com/brms/reference/brmsformula-helpers.html)
  (documented on
  [`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin_mfd_to_dpars.md)).
  The inverse reports `delta` as `NA` at `f = 1`, where there is no
  excess to be biased and the bias is genuinely unidentified, rather
  than silently returning `0`.

- New vignette, `Getting started with pairedcountbrms`: simulates paired
  counts from known `binegbin` parameters, fits them with
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) plus
  [`binegbin_stanvars()`](https://anhsmith.github.io/pairedcountbrms/reference/binegbin.md),
  checks that the five dpars recover the truth, and exercises
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html).
