# pairedcountbrms 0.7.0

* **Breaking (dpar rename): the joint families' excess-rate dpars
  `lambdaem`/`lambdalb` are now `lambdaone`/`lambdatwo`**, affecting `bipois()`,
  `binegbin()` and `binegbin_joint()`. The former names referred to electronic
  monitoring and vessel logbooks, the two sources in the project from which the
  package was extracted, and carried no meaning outside it. The new names are
  positional: source 1 is the response, source 2 is supplied via `vint()`.

  They are spelled out rather than written `lambda1`/`lambda2` because
  `brms::custom_family()` rejects dpar names ending in a digit, as well as dots
  and underscores. The documentation writes them as $\lambda_1$ and $\lambda_2$;
  a notation table in the README gives the correspondence.

* **Fits made before this release still post-process; no refitting is
  required.** A `brmsfit` stores its own family object, so
  `prepare_predictions()` on an older fit returns a `prep` carrying the previous
  dpar names. Rate reads in the joint families now resolve the name against the
  fit, so `log_lik()`, `loo()`, `posterior_predict()` and `posterior_epred()`
  return numerically identical results on pre-0.7.0 fits
  (`tests/testthat/test-dpar-compat.R`).

  Code that names the dpars does require editing: `nlf()`/`lf()` terms and
  `prior(..., dpar = )` arguments. Data column names are unaffected, since the
  second count and the observation flag reach the family positionally through
  `vint()`.

* **Breaking (argument names):** `binegbin_dpars_to_mfd()` takes
  `lambdaone`/`lambdatwo`, and `binegbin_mfd_to_dpars()` returns a list with
  those names. Positional calls are unaffected.

* **Documentation notation unified.** The joint families previously used
  $y_{\mathrm{em}}, y_{\mathrm{lb}}, N_{10}, N_{01}$ while the difference
  families used $y_1, y_2$, so the identity between the bivariate Poisson's
  difference and the Skellam distribution was stated in one notation and
  demonstrated in another. All documentation now uses $y_1, y_2$, $N_1, N_2$ and
  $\lambda_1, \lambda_2$, making $d = y_1 - y_2 = N_1 - N_2$ continuous across
  both suites.

* **Recovery tests separated into smoke gates and a calibration assessment.**
  The recovery tests asserted that the true value fell within a single fit's 90%
  credible interval. For a correct model with a calibrated posterior this is a
  Bernoulli(0.9) draw, failing 10% of the time by construction; across the ~18
  such assertions in the suite, spurious failures were the norm. Two assertions
  in `test-binegbin.R` had been failing on every run since they were written,
  identically before and after the dpar rename, undetected because
  `skip_on_cran()` excludes the fitting tests from `R CMD check`.

  Single-fit checks are now smoke gates: convergence plus a wide (99%) interval,
  which detects gross mis-specification without failing a correct model by
  chance. Calibration is assessed separately, as the proportion of repeated
  simulate-and-refit replicates whose nominal interval contains the truth, with
  the pass threshold derived from the Binomial(*R*, 0.9) null (0.16%
  false-failure rate at *R* = 10, against 10% per assertion previously).

  The coverage assessment costs minutes and is opt-in via
  `PAIREDCOUNTBRMS_COVERAGE=true`, a named gate rather than a silent default. It
  compiles the model once and reuses it through `update(recompile = FALSE)`,
  which keeps 10 replicates to approximately 3.5 minutes. Shared helpers are in
  `tests/testthat/helper-coverage.R`.

* **`binegbin` recovery and group-level composition are now separate tests, with
  less extreme generative truths.** A single test previously did both, on data
  with a vessel effect over 8 levels and severe overdispersion (the shared
  component's variance was 40 against a mean of 8). Three variance channels —
  the group-level SD, `shapes` and `shapex` — competed for the same residual.
  `shapex` is identified only through the difference's variance, and with the
  excess rate also unknown the two trade off along `shapex` = `lambda`^2 /
  (*V* − `lambda`) at fixed *V* = Var(*d*)/2, so a 1.5 SD fluctuation in Var(*d*)
  displaced the posterior substantially.

  Recovery is now asserted on an intercepts-only fit at *n* = 400 with
  approximately half of each component's variance from overdispersion
  (`shapes` = 8, `shapex` = 3), identifiable from both directions. A second test
  adds a group-level term and checks only that the family composes with brms's
  random-effects machinery and samples cleanly, making no dispersion-recovery
  claim. Recovery also includes a posterior-predictive check on Var(*d*), the
  quantity `shapex` governs.

* **References to the originating client project removed.** The documentation
  pointed readers at that project by name and at specific internal `.qmd` files
  — 26 references across the R sources, four help pages and three test files,
  several of which rendered onto the public documentation site, and none
  accessible to a reader. Where a pointer carried an argument (the
  regression-to-the-mean rationale for modelling the pair jointly, the
  observation-level random-effect failure behind the Negative-Binomial
  components, the prior-scale translation on `skellam1()`), that argument is now
  stated inline. Validation claims cite the package's own test files.

* **Continuous integration now runs the tests.** The repository previously
  contained one workflow, `pkgdown.yaml`, which built and deployed the
  documentation site and ran no tests. `R-CMD-check.yaml` adds two jobs on
  different cadences: `check` on every push and pull request, without a Stan
  backend (~5 minutes), and `check-stan` weekly and on demand, which sets
  `NOT_CRAN=true` and installs both rstan and cmdstanr so that the model fits
  execute.

  The fitting tests are both the most likely to detect a regression and too slow
  for every push. `skip_on_cran()` correctly excludes them from a routine check,
  but with nothing scheduled to lift it they had never run in continuous
  integration at all.

* **The $(M, f, \delta)$ reparameterisation is now demonstrated rather than
  described.** The anatomy article previously carried an unevaluated `bf()`
  fragment in place of a worked example — no `family`, no `stanvars`, no `brm()`
  call, and never executed. It now simulates from known coordinates, fits,
  recovers all five, and converts the posterior back through
  `binegbin_mfd_to_dpars()`.

  The fit is parameterised in $\kappa = 1/\sqrt{\phi}$ rather than $\phi$, via
  `nlf(shapes ~ -2 * log(kappas))`. In $\phi$ the Poisson limit is
  $\phi \to \infty$, so there is no finite point to shrink towards and an
  exponential prior on $\phi$ shrinks towards maximum overdispersion. In $\kappa$
  the limit is zero, which admits penalised-complexity priors (Simpson et al.
  2017): `exponential()` on each $\kappa$ shrinking to Poisson, and
  `double_exponential()` on $\delta$ shrinking to no between-source bias.
  Congruence $f$ is deliberately not shrunk in either direction, both of its
  ends being degenerate.

  Each prior is plotted on the scale it is interpreted on rather than the
  linear-predictor scale it is stated on.

* **The anatomy article is now precompiled**, as the vignette is, because it
  fits a model. `vignettes/articles/paired-count-anatomy.Rmd.orig` is the
  source; the `.Rmd` and `vignettes/articles/figure/` are generated. Edit the
  `.orig`. `vignettes/precompile.R` handles both documents, preserving the
  property that pkgdown requires no Stan backend.

* **New README section on priors.** `get_prior()` on a `binegbin()` model shows
  that brms assigns the `mu` dpar a default `student_t` prior but leaves `lamx`,
  `shapes` and `shapex` flat and improper; these are the weakly-identified
  parameters. In this package's recovery test, omitting priors produced a
  divergent transition and an Rhat of 1.0101. The section gives a
  weakly-informative set, notes that the `class`/`dpar`/`nlpar` slots differ
  between a rate supplied through `nlf()` and the dispersions, and advises
  shifting the prior mean rather than increasing its SD. The recovery tests use
  these priors.

* **New README note on `nl = TRUE`.** In a non-linear brms formula the main
  formula's right-hand side is a non-linear expression for `mu`, not a request
  for an intercept. `bf(y1 | vint(y2) ~ 1, nlf(...), ..., nl = TRUE)` without a
  `mu ~ 1` term therefore generates `mu[n] = exp(1)`, fixing the shared rate at
  *e*. This does not error; sampling becomes very slow. All the joint families'
  documented examples use `nl = TRUE` so that `nlf()` can tie the two excess
  rates.

* **`citation("pairedcountbrms")` now returns a usable citation.** There was no
  `inst/CITATION`, and the auto-generated fallback could not determine a year
  from `DESCRIPTION`, rendering as `Smith A (????)` with a warning.
  `inst/CITATION` reads the version and year from `DESCRIPTION`, and a `Date`
  field has been added. The header directs users to `citation("brms")` and to
  the distribution papers in the README.

* Examples and vignettes use `y1`/`y2` for the two count columns and `y1_obs`
  for `binegbin_joint()`'s observation flag, previously `y_em`/`y_lb`/`em_obs`.
  These are illustrative data column names, not API.

# pairedcountbrms 0.6.0

* **Renamed: `skellambrms` is now `pairedcountbrms`.** The old name named one
  family; the package's subject is the comparison of two paired count sources,
  by two complementary routes — the difference families (`skellam1`/`skellam2`,
  `dnorm1`/`dnorm2`, `dlaplace1`/`dlaplace2`) and the joint bivariate families
  (`bipois`, `binegbin`, `binegbin_joint`). Skellam is one member of that set,
  and no longer the most used one.

  **No family names change**, so no fitted model needs refitting: `binegbin`,
  `binegbin_joint`, `bipois`, `skellam1`/`skellam2`, `dnorm1`/`dnorm2` and
  `dlaplace1`/`dlaplace2` all keep their names, and brms continues to resolve
  each fit's `log_lik_*` / `posterior_predict_*` / `posterior_epred_*` methods
  off the attached search path exactly as before. The only change a user needs
  to make is `library(skellambrms)` → `library(pairedcountbrms)` (and any
  `skellambrms::` prefix).

  The GitHub repository moves to `anhsmith/pairedcountbrms`. GitHub serves a
  permanent redirect from the old path for both web and git, so existing
  clones and `pak::pak("anhsmith/skellambrms")` calls keep working.

* New coordinate helpers, `binegbin_mfd_to_dpars()` and
  `binegbin_dpars_to_mfd()`, converting between the rate dpars the joint
  families take (`mu`, `lambdaem`, `lambdalb`, plus `shapes`/`shapex`) and the
  interpretable `(M, f, delta)` coordinates — overall level, congruence, and
  method bias — along with the SD-scale dispersions `kappas`/`kappax`. Pure
  transforms; they fit nothing, and fitting in these coordinates still goes
  through `nlf()` (documented on `binegbin_mfd_to_dpars()`). The inverse
  reports `delta` as `NA` at `f = 1`, where there is no excess to be biased and
  the bias is genuinely unidentified, rather than silently returning `0`.

* New vignette, `Getting started with pairedcountbrms`: simulates paired counts
  from known `binegbin` parameters, fits them with `brm()` plus
  `binegbin_stanvars()`, checks that the five dpars recover the truth, and
  exercises `posterior_predict()` and `log_lik()`.

# skellambrms 0.5.0

* **Breaking (link change):** `binegbin()` and `binegbin_joint()` now log-link
  `lambdaem`/`lambdalb` (previously identity), so all five dpars are log-linked
  — the conventional log-linear rate parameterisation (Karlis & Ntzoufras 2003)
  and consistent with `bipois()`. brms applies the link *on top of* a
  non-linear formula, so the shared-excess idiom now drops the explicit
  `exp()`: write `nlf(lambdaem ~ lamx)` (the log link exponentiates) rather
  than `nlf(lambdaem ~ exp(lamx))`, which under the new link would
  double-exponentiate. The generated model, and hence the posterior, is
  identical to the old identity-link + explicit-`exp()` form; only the formula
  syntax changes. Migration: remove the `exp()` from every `nlf()` on
  `lambdaem`/`lambdalb`. Plain `lambdaem ~ 1` now gives a clean log-scale
  intercept instead of a bounded natural-scale one.
* Added `binegbin_joint()` / `binegbin_joint_stanvars()`: a **censoring-aware
  extension of `binegbin()`** for data in which the EM margin (`y_em`) is
  observed on only some rows. The same trivariate-reduction bivariate
  Negative-Binomial model and five dpars as `binegbin()`, but each row carries
  a second `vint()` integer, an `em_obs` 0/1 flag: on `em_obs == 1` (matched)
  rows the likelihood is the full joint `binegbin` lpmf on `(y_em, y_lb)`; on
  `em_obs == 0` (LB-only) rows it is the `y_em`-integrated marginal of that
  same joint (`sum_k NB2(k | mu, shapes) NB2(y_lb - k | lambdalb, shapex)`) --
  not a separate single-dispersion `neg_binomial_2` on `y_lb`. One `brm()`
  call thus pools matched and LB-only rows under one coherent likelihood:
  `lambdaem` and the EM/LB bias are identified only by the matched rows, while
  the LB-only rows sharpen `mu`, `shapes`, `lambdalb`, and the shared
  vessel/trip random-effect structure.
* The second count and the flag travel via `vint(y_lb, em_obs)`
  (`vint1 = y_lb`, `vint2 = em_obs`). `binegbin_joint()` exports the standard
  `log_lik_binegbin_joint()` and `posterior_predict_binegbin_joint()`
  interface functions; the latter simulates `y_em` conditional on the observed
  `y_lb` for **every** row (matched and LB-only alike, ignoring `em_obs`) via
  the same discrete `N_shared | y_lb` conditional as `binegbin()`.
* On `em_obs == 1` rows `binegbin_joint`'s lpmf equals `binegbin`'s exactly;
  the suite pins this equivalence (R reference and Stan), the marginal
  identity (sum over `y_em` of the matched branch == the LB-only branch), and
  the conditional-prediction identity (`posterior_predict` draws ==
  joint / marginal), alongside normalisation, a Stan-vs-R grid cross-check to
  ~1e-14, numerical-stability edge cases, and a censored brms end-to-end fit
  that verifies `loo()`/`posterior_predict()` dispatch and parameter recovery.
* Truncation (`resp_trunc()`) is not applicable to this joint family and no
  `_lccdf_stanvars()` is provided.

# skellambrms 0.4.0

* Added the first **joint bivariate-count families**, a different modelling
  paradigm from the six difference families above: instead of modelling
  `y_em - y_lb` (a Z-valued difference), these model the matched pair
  `(y_em, y_lb)` jointly, capturing their correlation and marginal
  overdispersion alongside their difference. Both are built by trivariate
  reduction — `y_em = N_shared + N10`, `y_lb = N_shared + N01`, with
  `N_shared` marginalised out of the joint likelihood analytically — and take
  the second count via brms's `vint()` addition term.
* Added `bipois()` / `bipois_stanvars()`: the **bivariate Poisson**, with
  three independent Poisson latent components (`mu` = shared rate,
  `lambdaem`/`lambdalb` = the two private rates). The Stan log-likelihood uses
  the incremental `bipois2` recurrence (stan-users, March 2016), cross-checked
  against an independent R brute-force reference to ~1e-14. Cannot be
  overdispersed (each component has `Var == mean`), so it underfits any data
  whose margins are overdispersed — see `binegbin()`.
* Added `binegbin()` / `binegbin_stanvars()`: the **bivariate
  Negative-Binomial**, the overdispersed sibling of `bipois()`. Each latent
  component is Negative-Binomial (`neg_binomial_2`), adding two scalar
  dispersion dpars — `shapes` (shared component) and `shapex` (shared across
  the two private components, the "2 kappa" structure). The marginalisation
  sum is identical in form to `bipois()` with `neg_binomial_2_lpmf` swapped in
  for `poisson_lpmf` (not a Gamma-mixed Poisson, so no stacked
  marginalisation). Carries overdispersion in identifiable *scalar*
  parameters rather than per-observation random effects, which was found to
  overfit (the excess-dispersion SD collapses under a per-set OLRE). Validated
  against an independent R brute-force reference to ~1e-14, with normalisation
  and moment-identity checks and clean synthetic parameter recovery.
* Both joint families export the standard `log_lik_<family>()`,
  `posterior_predict_<family>()`, and `posterior_epred_<family>()` interface
  functions. `posterior_predict` simulates `y_em` conditional on the observed
  `y_lb` — for `bipois()` via the closed-form `Binomial(y_lb, mu/(mu+lambdalb))`
  split, for `binegbin()` via the discrete `N_shared | y_lb` conditional
  (no clean Binomial form for a NegBin sum). Truncation (`resp_trunc()`) is
  not applicable to these joint families and no `_lccdf_stanvars()` is
  provided.

# skellambrms 0.3.2

* Fixed a silent `ifelse()` length-collapse bug in `log_lik_dlaplace1()`:
  `ifelse(test, yes, no)` takes its output length from `test`, not from
  the (vectorised) `yes`/`no` branches. Because `dlaplace1()` has no free
  mean, the CDF-differencing argument built from the observation `z` is a
  scalar, while `sigma` (and the derived `b`) varies across posterior
  draws — so the `ifelse()` test was evaluated at length 1 and the whole
  per-observation log-likelihood silently collapsed to length 1 instead of
  `ndraws`. This broke `brms::add_criterion(fit, "loo")` for every
  `dlaplace1()` fit (`is.matrix(unnormalized_log_weights) is not TRUE`),
  while sampling, `posterior_predict()`, and truncation were entirely
  unaffected — confirmed isolated to this one function's R-side length
  handling, not a data or convergence issue.
* Applied the same fix pre-emptively to five more internal R-side helpers
  sharing the identical `ifelse()` shape — `dlaplace1_lpmf_r()`,
  `dlaplace1_lccdf_r()`, `dlaplace2_lccdf_r()`, `dnorm1_lpmf_r()`, and
  `dnorm2_lpmf_r()` in `R/truncation.R`. None were triggering the bug at
  their current call sites (which happen to keep argument lengths
  matched), but all shared the same landmine.

# skellambrms 0.3.1

* Fixed `posterior_predict_<family>()` for all six families: previously
  ignored `trunc()`/`resp_trunc()` bounds entirely, drawing from the
  untruncated distribution and returning it verbatim (confirmed to produce
  out-of-bound draws, e.g. `posterior_predict_dnorm2()` returning values
  well below a `lb = -14` bound). Now performs correct inverse-CDF sampling
  within the truncation bounds, reusing each family's already-validated
  log-CCDF math (transcribed to R in the new internal `R/truncation.R`)
  rather than rejection sampling, which was confirmed empirically
  slow/low-acceptance for tight bounds — especially costly for
  `skellam2()`, whose per-evaluation cost (an iterative Bessel-function
  tail-sum) is comparatively high.
* Fixed `posterior_epred_<family>()` for all six families: previously
  returned the untruncated mean (`mu`, or `0` for `skellam1()`/
  `dlaplace1()`/`dnorm1()`) even when a truncation bound was tight enough
  to meaningfully shift the conditional expectation, with no warning. Now
  computes the correct truncated conditional expectation via deterministic
  numerical summation of the truncated PMF — not Monte Carlo, so the
  result is exact to a documented tolerance and fully reproducible.
* **Known limitation surfaced (not introduced) by this fix:**
  `brms::posterior_epred()` — and anything built on it, including
  `fitted()` and `conditional_effects()` — errors on any truncated fit of
  a custom family, for all six families here, under the currently
  installed `brms`. This is a `brms` limitation: its internal dispatcher
  checks whether a fit is truncated *before* checking whether the family
  is a custom one, and has no fallback to a custom family's own
  `posterior_epred_<family>()` on the truncated branch. Call
  `posterior_epred_<family>(brms::prepare_predictions(fit))` directly as a
  workaround. `brms::posterior_predict()` is unaffected and works
  correctly for truncated fits of every family. See the README's
  "Limitations" section for details.

# skellambrms 0.3.0

* **Breaking change:** `skellam1()` now samples on `sigma`, the SD of the
  difference (log-linked), rather than the underlying Skellam rate
  directly. `mu_skellam = sigma^2 / 2` is derived internally; the
  Bessel-sum likelihood itself is unchanged. A prior previously stated on
  `log(mu_skellam)` translates as
  `log(sigma) = 0.5*log(2) + 0.5*log(mu_skellam)` — e.g. an old
  `normal(1, 1.5)` becomes `normal(0.847, 0.75)`. This reparameterisation
  establishes a common (mean, SD-scale) convention shared by every family
  below.
* Added `skellam2()` / `skellam2_stanvars()` / `skellam2_lccdf_stanvars()`:
  the asymmetric Skellam (Koopman parameterisation), with a free mean
  (`mu`) and a free `sigmaexcess` (so that
  `sigma^2 = |mu| + sigmaexcess^2`, guaranteeing Skellam validity for
  every `mu` and `sigmaexcess >= 0` — a corrected constraint relative to
  the originally-specified `sigma = sqrt(mu^2 + sigmaexcess^2)`, which
  only guarantees the weaker `sigma >= |mu|` and admits invalid
  (negative-rate) parameter combinations for `|mu| < 1`). Reduces exactly
  to `skellam1()` at `mu = 0`.
* Added `dlaplace1()` / `dlaplace1_stanvars()` / `dlaplace1_lccdf_stanvars()`:
  a discrete Laplace distribution (location fixed at 0, free `sigma`),
  discretised from the continuous Laplace via CDF differencing.
* Added `dlaplace2()` / `dlaplace2_stanvars()` / `dlaplace2_lccdf_stanvars()`:
  the free-location/free-scale discrete Laplace, with no constraint
  coupling `mu` and `sigma` — a deliberate structural contrast with
  `skellam2()`, for comparing models where bias and spread are
  structurally coupled against ones where they vary independently.
* Added `dnorm1()` / `dnorm1_stanvars()` / `dnorm1_lccdf_stanvars()`: a
  discrete normal distribution (location fixed at 0, free `sigma`), via
  the same CDF-differencing pattern as `dlaplace1()`.
* Added `dnorm2()` / `dnorm2_stanvars()` / `dnorm2_lccdf_stanvars()`: the
  free-location/free-scale discrete normal, structurally analogous to
  `dlaplace2()`.
* Fixed a numerical-stability issue affecting `skellam1_lccdf_stanvars()`
  and `skellam2_lccdf_stanvars()`'s normal-approximation branch, and
  `dnorm1`/`dnorm2`'s `_lpmf`/`_lccdf`: Stan's built-in `normal_lccdf` is
  not safe to call directly in this context. This is a documented Stan
  limitation, not a guess — the Stan Functions Reference states
  `normal_lccdf` underflows to `-inf` for `(y-mu)/sigma > ~8.25`, and
  [stan-dev/math#1985](https://github.com/stan-dev/math/issues/1985)
  confirms `normal_lccdf` (unlike `normal_lcdf`) was never updated with
  the more accurate Mills-ratio approximation. Fixed via an exact
  `erfc()`-based closed form throughout, confirmed to match a trusted R
  reference to machine precision out to 30+ SDs.
* Fixed a Stan-compiler portability bug: `skellam2_lpmf`/`skellam2_lccdf`
  used `fabs()`, which compiles under `rstan`'s bundled Stan version but
  is not a valid identifier under `cmdstanr`'s (use `abs()`, which is
  type-generic and already used elsewhere in the same functions).
* Added `cmdstanr` to `Suggests` (previously only `rstan` was declared,
  so `R CMD check`'s isolated test environment could not see an
  already-installed `cmdstanr`).

# skellambrms 0.2.0

* Added `skellam1_lccdf_stanvars()`, providing the log-CCDF of the
  symmetric Skellam(mu, mu) distribution so that brms's `resp_trunc()`
  can be used with `skellam1()`, including row-varying truncation
  bounds. Still the symmetric Skellam(mu, mu) case only — this adds
  truncation support to the existing family, not a new family or the
  asymmetric case.
* The exact log-CCDF (an iterative Bessel-sum tail) switches to a normal
  approximation above a configurable `normal_approx_threshold` (default
  `100`), guarding against a confirmed `std::bad_alloc` crash and a
  confirmed multi-GB memory blowup when an unadapted HMC proposal pushes
  `mu` to an extreme value during warmup. See `?skellam1_lccdf_stanvars`
  for guidance on choosing this threshold for your own data.

# skellambrms 0.1.0

* Initial release: `skellam1()` and `skellam1_stanvars()`, a brms custom
  family for the symmetric Skellam(mu, mu) distribution.
