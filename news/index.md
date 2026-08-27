# Changelog

## bicountbrms 0.10.0

- **Breaking: four families become two, with two constructors each.**
  The package now supplies
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  and
  [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md),
  each with a `_partialobs()` sibling:

  |  | fully paired | first count missing on some rows |
  |----|----|----|
  | Poisson | [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md) | [`bipois_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/bipois_partialobs.md) |
  | Negative-Binomial | [`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md) | [`binegbin_partialobs()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_partialobs.md) |

  Both constructors of a pair return the *same* `custom_family` name, so
  there is one Stan `_lpmf` and one set of `log_lik_*` /
  `posterior_predict_*` / `posterior_epred_*` methods per component
  distribution. They differ only in whether the observation flag reaches
  the likelihood as data or as a constant. A user fitting fully paired
  data never has to supply the flag, or know it exists.

  `bipois_cens()`, `binegbin_cens()` and the five `binegbin_joint` names
  are **removed outright**. There is no deprecation layer and nothing
  left to remove at a later version; see the last bullet for what to do
  if you hold a fit made under those names.

- **[`binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  now has six dpars.** `shapex` is gone; both excess components have
  their own dispersion, `shapexone` and `shapextwo`, as
  `binegbin_cens()` has had since 0.8.0. This is the change the release
  exists for. The capability was on the wrong family: in a partially
  paired design the unmatched rows inform only `shapextwo`, so
  `shapexone` is identified by the matched rows alone, whereas in a
  fully paired design every row informs both and they are as easy to
  identify as they ever get.

  The symmetric model is a formula constraint, exactly as it has been
  for the partially observed family since 0.8.0:

  ``` r

  bf(y1 | vint(y2) ~ 1,
     nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
     shapexx ~ 1, ..., nl = TRUE)
  ```

  Set priors on it with `nlpar = "shapexx"` rather than
  `dpar = "shapex"`.

- **`_cens` named the wrong mechanism.** Censoring means a value is
  known to lie in a set. On these rows the first count is not observed
  at all and the likelihood marginalises over its whole support. brms
  already uses `cens()` as an addition term meaning exactly the bounded
  thing — `left`, `right`, `interval`, with no code for “not observed” —
  so a brms user met a familiar word attached to an incompatible
  mechanism. `partialobs` names what is actually partial: the *pair*,
  which is a property of the data. `unobs` and `unmatched` name only the
  minority row class, and `partial` alone collides with partial
  likelihood and partial pooling.

- **One Stan function for both constructors, through a literal in
  `vars`.** brms builds the generated lpmf call by pasting each entry of
  `family$vars` in verbatim, and validates those entries no further than
  [`as.character()`](https://rdrr.io/r/base/character.html). The plain
  constructors exploit that: they declare `vars = c("vint1[n]", "1")`,
  so a fully paired model calls the same function with `y1_obs` fixed at
  the literal `1`.

      binegbin()             target += binegbin_lpmf(Y[n] | ..., vint1[n], 1);
      binegbin_partialobs()  target += binegbin_lpmf(Y[n] | ..., vint1[n], vint2[n]);

  The alternative — two Stan functions of the same name and different
  arity — needs user-defined function overloading, which arrived in Stan
  2.29 (February 2022). Taking it would have obliged `DESCRIPTION` to
  declare floors on rstan and cmdstanr, whose absence would otherwise
  surface as an opaque compile error. This design needs none, and none
  has been added.

  That pasting behaviour is not a documented brms guarantee, so
  `tests/testthat/test-stancode-shape.R` pins it: it asserts the exact
  call each constructor generates and that the one-`vint`
  [`standata()`](https://paulbuerkner.com/brms/reference/standata.html)
  contains no `vint2`. It needs brms but no Stan toolchain, so it runs
  in the fast suite.

- **New tests for the fully paired path.** Most users take the fully
  paired path and, before this release, nothing tested its prediction,
  expectation or dispersion routing: with a single `shapex` there was
  nothing to get wrong.

  - `tests/testthat/test-partialobs-predict.R` replaces
    `test-cens-predict.R`. It keeps all five of that file’s blocks and
    adds a one-`vint` mirror of each: the same fixtures with `shapexone`
    and `shapextwo` an order of magnitude apart, on a prep with no
    `vint2`.
  - `tests/testthat/test-unified-vint.R` pins the release criterion that
    the two shapes agree: for each of `log_lik`, `posterior_predict` and
    `posterior_epred`, the one-`vint` path on matched data equals the
    two-`vint` path with the flag set to 1. `log_lik` is the only method
    that reads the flag, so its non-vacuity partner asserts the flag-0
    answer *differs*; the other two ignore the flag by design, and that
    is asserted rather than left implicit.
  - `tests/testthat/test-dpar-compat.R` gains the stored-fit shape
    (below).

- **Stored `binegbin` fits keep working.** A fit stores its own family
  object, and brms resolves post-processing from the name it stored, so
  a five-dpar `binegbin` fit lands on the new six-dpar methods however
  this package is pinned. Three independent fallbacks have to line up
  for that to work: `.get_rate()` resolves pre-0.7.0
  `lambdaem`/`lambdalb`; `.SHAPEXONE_NAMES` and `.SHAPEXTWO_NAMES` both
  list `shapex` last, so both per-margin dispersions resolve to the one
  such a fit has; and an absent `vint2` selects the matched branch.
  Scanning the sibling project’s fit directory found 103 stored
  `binegbin` fits, every one of them with the pre-0.7.0 rate names *and*
  the single `shapex`, so all three paths are exercised in practice.
  `test-dpar-compat.R` now builds that exact prep and asserts all three
  methods agree with the six-dpar answer with the dispersions tied.

- **[`binegbin_mfd_to_dpars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_mfd_to_dpars.md)
  no longer emits `shapex`.** Supplying `kappax` now writes `shapexone`
  and `shapextwo` at that common value — the tied, symmetric model —
  rather than a single `shapex`, which is no longer a dpar of anything
  and so could not be passed to
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html). `kappax`
  survives as the shorthand it always was; only the output naming
  changed.

  [`binegbin_dpars_to_mfd()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_dpars_to_mfd.md)
  still *accepts* `shapex`, because its input is a stored fit and that
  is the genuine dpar name on any `binegbin` fit made before this
  release. The two directions differ deliberately: the forward one
  writes the current vocabulary, the inverse reads the old one, which is
  the same principle `.get_dpar_any()` applies.

- **The README is now an orientation page.** It ran to about 5,700
  words, which is not a length anyone reads, and it set out the
  construction, the moments, prior guidance, validation and the release
  history. It is now about 800: what the package is, how to install it,
  the four constructors in a table, the notation bridge from symbol to
  distributional-parameter name, one minimal fit, and links.

  Four articles are new, and the site gains an `articles:` section,
  which it did not have before: neither the get-started vignette nor
  `paired-count-anatomy` was previously placed on it deliberately.

  - *The families and their parameters* — the construction, the moments,
    the native distributional parameters, and a fitted demonstration
    that two excess dispersions an order of magnitude apart are resolved
    as different rather than merely bracketed.
  - *Choosing priors* — the improper-defaults fact, the recipes, and
    which of the `class`, `dpar` and `nlpar` slots each prior belongs
    in. It routes on to `paired-count-anatomy` for shrinkage priors in
    the interpretable coordinates.
  - *A worked partially observed fit* — the first such example in the
    package. A design with 40% of first counts unrecorded is simulated,
    fitted, and the imputed first count scored against the values
    withheld from the model.
  - *Migration and errata* — what each release changed, the
    substitutions for existing code, and the two corrected numerical
    results.

- **The recommended priors in the $`(M, f, \delta)`$ coordinates have
  changed.** *The anatomy of a paired count* now recommends
  `normal(0, 0.5)` on `methd`, and `normal(0, 1)` with `lb = 0`, a
  half-normal, on each of `kappas` and `kappax`. The previous
  recommendation was `double_exponential(0, 0.5)` and `exponential(1)`.
  The worked fit in that article was refitted under the new priors:
  every 90% interval still covers its true value, and no posterior
  median moved by more than 1%. `exponential(1)` is documented as an
  alternative of the same shape with a heavier tail.

- **Those priors are no longer called penalised-complexity priors.**
  They follow the same principles: a base model at zero, and a density
  decreasing monotonically away from it. The PC construction
  additionally requires a distance measure derived for the specific
  model and a constant-rate assumption on that distance, neither of
  which is established for these families. The wording is corrected in
  `README.md`, the get-started vignette, *Choosing priors* and *The
  anatomy of a paired count*.

- **If you hold a fit made with `bipois_cens()`, `binegbin_cens()` or
  `binegbin_joint()`**, install the last version that defined them —
  0.9.1 — and keep it on the search path for that fit. brms resolves
  post-processing late, by looking up `log_lik_<name>` from the stored
  family name at call time, so
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  need the methods of the version the fit was made under. No argument to
  those calls can route around it, and 0.10.0 deliberately ships no
  shim: pinning the whole package is a stronger guarantee than a
  hand-picked subset of forwarders, and the one project with such fits
  already pins 0.9.1 in its own `renv.lock`. \# bicountbrms 0.9.1

- **No behaviour change.** Nothing in the four families’ likelihoods,
  predictions or expectations differs from 0.9.0. This release adds test
  coverage that was missing and corrects a claim in the README the suite
  did not support. No fitted model is affected and no code needs
  changing.

- **New tests: prediction on censored rows, with the two excess
  dispersions distinct** (`tests/testthat/test-cens-predict.R`). The
  `_cens` families exist to impute the first margin on rows where it was
  never observed, and that row class had no correctness test. The Monte
  Carlo check validating `posterior_predict` against the exact
  conditional `P(y1 | y2)` ran at `y1_obs == 1` in both
  `test-binegbin_cens.R` and `test-bipois_cens.R` — matched rows, the
  one class where it cannot fail. Every prediction test also built its
  prep from a single five-dpar `shapex`, which `.SHAPEXONE_NAMES` and
  `.SHAPEXTWO_NAMES` both resolve to, so the two dispersions were never
  distinguished from each other.

  That combination left a specific defect invisible. In
  `posterior_predict_binegbin_cens()` the conditional split of
  `N_shared | y2` is weighted by `shapextwo` while the fresh `N1` added
  on top takes `shapexone`; exchanging them passed 412 assertions across
  six test files with no failures. The new file fails on it four times.
  It also pins that `y1_obs` does not change the draws — the flag
  selects a likelihood branch, not a prediction — and that `y2 = 0`
  leaves `posterior_predict` drawing `NB2(lambdaone, shapexone)` alone,
  which is the short-circuit keeping
  [`sample()`](https://rdrr.io/r/base/sample.html) from being handed an
  empty support.

- **README: the conditional-prediction identity is now claimed
  accurately.** The Testing section asserted it for the `_cens` families
  without qualification while the suite checked it on matched rows only.
  It now states that it holds on matched *and* censored rows and with
  the two excess dispersions distinct, which is what the tests above
  pin.

## bicountbrms 0.9.0

- **The package has been split, and renamed.** `pairedcountbrms` carried
  two unrelated suites: families that model a count pair jointly, and
  families that model its difference. They shared no code — not a
  helper, not a Stan block, not a dependency. This release separates
  them. The joint families are `bicountbrms`; the difference families
  (`skellam1`/`skellam2`, `dnorm1`/`dnorm2`, `dlaplace1`/`dlaplace2`)
  have returned to
  [`skellambrms`](https://github.com/anhsmith/skellambrms), the name
  under which they were originally released.

  **No fitted model needs refitting**, either for the split or for the
  rename below. brms resolves each fit’s `log_lik_*` /
  `posterior_predict_*` / `posterior_epred_*` methods off the attached
  search path at call time, so a stored fit works as soon as the package
  holding its family is attached. No dpar name changed. The only source
  change the split itself requires is
  [`library(pairedcountbrms)`](https://rdrr.io/r/base/library.html) →
  [`library(bicountbrms)`](https://github.com/anhsmith/bicountbrms) or
  [`library(skellambrms)`](https://rdrr.io/r/base/library.html), and any
  `pairedcountbrms::` prefix.

  The GitHub repository is renamed to `anhsmith/bicountbrms`, and GitHub
  serves a permanent redirect from `anhsmith/pairedcountbrms` for both
  web and git, so existing clones and
  `pak::pak("anhsmith/pairedcountbrms")` keep working. The redirect from
  `anhsmith/skellambrms` no longer applies: that name is now an
  independent repository holding the difference families, which is what
  a user who installed `skellambrms` was asking for.

  The timing is deliberate. A forthcoming 1.0.0 will be archived and
  assigned a DOI, and cited from a methods article describing the joint
  families; the citable artefact should therefore describe one thing
  rather than half of two. This 0.9.0 makes the split and the additions
  below available for checking first, so that what is archived has been
  exercised rather than merely tested.

- **Renamed: `binegbin_joint()` → `binegbin_cens()`.** All four families
  in this package model the two counts jointly, so `_joint` named a
  property common to all of them. The property that selects this one is
  that the first margin may be unobserved on some rows, which `_cens`
  names. `binegbin_cens_stanvars()`, `binegbin_cens_lpmf` (Stan) and the
  three post-processing methods follow the same rename.

  **Old code and old fits both keep working.** `binegbin_joint()` and
  `binegbin_joint_stanvars()` remain as deprecated wrappers that warn
  once per call and return the new objects. The three post-processing
  functions `log_lik_binegbin_joint()`,
  `posterior_predict_binegbin_joint()` and
  `posterior_epred_binegbin_joint()` remain as *silent* forwarders,
  because a fit made before this release stores
  `family$name == "binegbin_joint"` permanently and brms builds its
  method names from that stored name — nothing the caller passes to
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) or
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  can route around it. All five are removed in the next major version.

  No dpar name changed, so the forwarding is a straight hand-off, and a
  five-dpar pre-0.8.0 fit still resolves through both compatibility
  layers at once (a test pins this).

- **New family: `bipois_cens()` / `bipois_cens_stanvars()`.** The
  censoring-aware bivariate Poisson — `binegbin_cens()`’s equidispersed
  counterpart, and
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)’s
  censoring-aware one. Three log-linked dpars (`mu`, `lambdaone`,
  `lambdatwo`) and two `vint()` integers, `vint(y2, y1_obs)`, exactly as
  `binegbin_cens()` takes them. It has no deprecated aliases: it was
  added in this release under its final name.

  Its censored branch is closed form. A sum of independent Poissons is
  Poisson, so the `y1`-integrated marginal collapses to
  `poisson_lpmf(y2 | mu + lambdatwo)` — no convolution sum, no cutoff.
  `binegbin_cens()` cannot do this (NB2 + NB2 is not NB2 unless the two
  components share a success probability), which is why it evaluates
  that branch term by term.

  Two things follow. Users whose counts are equidispersed no longer have
  to fit `binegbin_cens()` with its dispersions pressed against the
  Poisson boundary, where sampling degrades. And the suite gains an
  independent analytic reference for `binegbin_cens()`’s censored
  branch: driving its three dispersions to their Poisson limit must
  reproduce `bipois_cens()`, checked as a limit — the error must shrink
  in proportion to `1/phi` — rather than at a single tolerance.
  Previously that branch was pinned only by the marginal identity, which
  compares it against the matched branch of the same code and so shares
  any error common to both sums.

- **Breaking (numerical):
  [`posterior_epred_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  is now exact.** It previously substituted the *marginal* shared
  fraction `mu/(mu + lambdatwo)` for the *conditional* one — the
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)
  answer, which is exact only in the Poisson limit — and its docstring
  described this as a point approximation adequate for display. It is
  not adequate. Over a grid of plausible rates and dispersions, with
  `y2` within one standard deviation of its mean, the substitution was
  in error by more than 5% in about half the settings and more than 20%
  in a third, in either direction depending on which component carried
  the greater dispersion.

  `E[N_shared | y2]` is now evaluated as the mean of the discrete
  conditional over `k = 0..y2` — the same weights
  [`posterior_predict_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  already samples from, summed rather than sampled — so the expectation
  and the predictions agree by construction rather than approximately.

  Numbers produced by an earlier version’s
  [`posterior_epred_binegbin()`](https://anhsmith.github.io/bicountbrms/reference/binegbin.md)
  should be recomputed.
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  are unaffected; the likelihood has not changed.

- **New: `posterior_epred_binegbin_cens()`.** That family previously
  defined no `posterior_epred` at all, on the view that `E[y1]` is
  ambiguous when `y1` is unobserved. It is not: `E[y1 | y2]` is well
  defined whether or not `y1` was recorded, and it is the quantity a
  user imputing the missing margin wants. It is returned for every row,
  matched and censored alike, matching
  `posterior_predict_binegbin_cens()`’s existing convention. Stored fits
  reach it through the deprecated `posterior_epred_binegbin_joint()`
  forwarder, so they gain the method without refitting.

- **All four families now share one `posterior_epred` convention**,
  stated once in `tests/testthat/test-epred.R` and held to one standard:
  the expectation must equal the mean of that family’s own
  `posterior_predict` draws to within Monte Carlo error. Before this
  release the three families answered three different ways — exactly,
  approximately, and not at all.

- `skellam` moves from **Imports to Suggests**. No function in `R/`
  calls it now that the Skellam families have left;
  `tests/testthat/test-bipois.R` still uses
  [`skellam::dskellam()`](https://rdrr.io/pkg/skellam/man/skellam.html)
  as the independent reference for the induced difference distribution,
  which is worth keeping.

- The coverage gate’s environment variable is renamed
  `PAIREDCOUNTBRMS_COVERAGE` → `BICOUNTBRMS_COVERAGE`.

- The README and the getting-started vignette are rewritten around the
  four joint families, and the vignette’s “one limitation to know about”
  section — which described brms’s
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  failure on *truncated* custom-family fits — is replaced by a statement
  of the epred convention above. That limitation cannot arise here:
  these families do not support
  [`resp_trunc()`](https://paulbuerkner.com/brms/reference/addition-terms.html).
  It applies to the difference families and is now documented in
  `skellambrms`.

------------------------------------------------------------------------

Entries below record this code’s history under its former names. It was
`skellambrms` through 0.5.0 and `pairedcountbrms` from 0.6.0 to 0.8.0,
and those packages contained the difference families as well, so some
entries describe families that are no longer here. They are left as
written rather than retrospectively edited. Releases before 0.4.0, which
predate the joint families entirely, are recorded in [`skellambrms`’s
NEWS](https://github.com/anhsmith/skellambrms/blob/master/NEWS.md).
