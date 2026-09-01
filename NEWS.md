# bicountbrms 0.10.0

* **Breaking: four families become two, with two constructors each.** The
  package now supplies `bipois()` and `binegbin()`, each with a
  `_partialobs()` sibling:

  |  | fully paired | first count missing on some rows |
  |---|---|---|
  | Poisson | `bipois()` | `bipois_partialobs()` |
  | Negative-Binomial | `binegbin()` | `binegbin_partialobs()` |

  Both constructors of a pair return the *same* `custom_family` name, so there
  is one Stan `_lpmf` and one set of `log_lik_*` / `posterior_predict_*` /
  `posterior_epred_*` methods per component distribution. They differ only in
  whether the observation flag reaches the likelihood as data or as a constant.
  A user fitting fully paired data never has to supply the flag, or know it
  exists.

  `bipois_cens()`, `binegbin_cens()` and the five `binegbin_joint` names are
  **removed outright**. There is no deprecation layer and nothing left to
  remove at a later version; see the last bullet for what to do if you hold a
  fit made under those names.

* **`binegbin()` now has six dpars.** `shapex` is gone; both
  excess components have their own dispersion, `shapexone` and `shapextwo`, as
  `binegbin_cens()` has had since 0.8.0. This is the change the release exists
  for. The capability was on the wrong family: in a partially paired design
  the unmatched rows inform only `shapextwo`, so `shapexone` is identified by the
  matched rows alone, whereas in a fully paired design every row informs both
  and they are as easy to identify as they ever get.

  The symmetric model is a formula constraint, exactly as it has been for the
  partially observed family since 0.8.0:

  ```r
  bf(y1 | vint(y2) ~ 1,
     nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx),
     shapexx ~ 1, ..., nl = TRUE)
  ```

  Set priors on it with `nlpar = "shapexx"` rather than `dpar = "shapex"`.

* **`_cens` named the wrong mechanism.** Censoring
  means a value is known to lie in a set. On these rows, the first count is not
  observed at all and the likelihood marginalises over its whole support. brms
  already uses `cens()` as an addition term meaning exactly the bounded thing —
  `left`, `right`, `interval`, with no code for "not observed" — so a brms user
  met a familiar word attached to an incompatible mechanism. `partialobs` names
  what is actually partial: the *pair*, which is a property of the data.
  `unobs` and `unmatched` name only the minority row class, and `partial` alone
  collides with partial likelihood and partial pooling.

* **One Stan function for both constructors, through a literal in `vars`.**
  brms builds the generated lpmf call by pasting each entry of
  `family$vars` in verbatim, and validates those entries no further than
  `as.character()`. The plain constructors exploit that: they declare
  `vars = c("vint1[n]", "1")`, so a fully paired model calls the same function
  with `y1_obs` fixed at the literal `1`.

  ```
  binegbin()             target += binegbin_lpmf(Y[n] | ..., vint1[n], 1);
  binegbin_partialobs()  target += binegbin_lpmf(Y[n] | ..., vint1[n], vint2[n]);
  ```

  The alternative — two Stan functions of the same name and different arity —
  needs user-defined function overloading, which arrived in Stan 2.29
  (February 2022). Taking it would have obliged `DESCRIPTION` to declare floors
  on rstan and cmdstanr, whose absence would otherwise surface as an opaque
  compile error. This design needs none, and none has been added.

  That pasting behaviour is not a documented brms guarantee, so
  `tests/testthat/test-stancode-shape.R` pins it: it asserts the exact call
  each constructor generates and that the one-`vint` `standata()` contains no
  `vint2`. It needs brms but no Stan toolchain, so it runs in the fast suite.

* **New tests for the fully paired path.** Most users take the fully paired
  path and, before this release, nothing tested its
  prediction, expectation or dispersion routing: with a single `shapex` there
  was nothing to get wrong.

  - `tests/testthat/test-partialobs-predict.R` replaces `test-cens-predict.R`.
    It keeps all five of that file's blocks and adds a one-`vint` mirror of
    each: the same fixtures with `shapexone` and `shapextwo` an order of
    magnitude apart, on a prep with no `vint2`.
  - `tests/testthat/test-unified-vint.R` pins the release criterion that the
    two shapes agree: for each of `log_lik`, `posterior_predict` and
    `posterior_epred`, the one-`vint` path on matched data equals the
    two-`vint` path with the flag set to 1. `log_lik` is the only method that
    reads the flag, so its non-vacuity partner asserts the flag-0 answer
    *differs*; the other two ignore the flag by design, and that is asserted
    rather than left implicit.
  - `tests/testthat/test-dpar-compat.R` gains the stored-fit shape (below).

* **Stored `binegbin` fits keep working.** A fit stores its own family
  object, and brms resolves
  post-processing from the name it stored, so a five-dpar `binegbin` fit lands
  on the new six-dpar methods however this package is pinned. Three independent
  fallbacks have to line up for that to work: `.get_rate()` resolves pre-0.7.0
  `lambdaem`/`lambdalb`; `.SHAPEXONE_NAMES` and `.SHAPEXTWO_NAMES` both list
  `shapex` last, so both per-margin dispersions resolve to the one such a fit
  has; and an absent `vint2` selects the matched branch. Scanning the sibling
  project's fit directory found 103 stored `binegbin` fits, every one of them
  with the pre-0.7.0 rate names *and* the single `shapex`, so all three
  paths are exercised in practice. `test-dpar-compat.R` now builds that exact
  prep and asserts all three methods agree with the six-dpar answer with the
  dispersions tied.

* **`binegbin_mfd_to_dpars()` no longer emits `shapex`.** Supplying `kappax`
  now writes `shapexone` and `shapextwo` at that common value — the tied,
  symmetric model — rather than a single `shapex`, which is no longer a dpar of
  anything and so could not be passed to `brm()`. `kappax` survives as the
  shorthand it always was; only the output naming changed.

  `binegbin_dpars_to_mfd()` still *accepts* `shapex`, because its input is a
  stored fit and that is the genuine dpar name on any `binegbin` fit made
  before this release. The two directions differ deliberately: the forward one
  writes the current vocabulary, the inverse reads the old one, which is the
  same principle `.get_dpar_any()` applies.

* **The README is now an orientation page.** The construction, the moments, prior
  guidance, validation and the release history have moved into the articles. What
  remains is what the package is, how to install it, the constructors in a table,
  the notation bridge from symbol to distributional-parameter name, one minimal
  fit, and links.

  Four articles are new, and the site gains an `articles:` section, which it
  did not have before: neither the get-started vignette nor
  `paired-count-anatomy` was previously placed on it deliberately.

  - *The families and their parameters* — the construction, the moments, the
    native distributional parameters, and a fitted demonstration that two
    excess dispersions an order of magnitude apart are resolved as different
    rather than merely bracketed.
  - *Choosing priors* — the improper-defaults fact, the recipes, and which of
    the `class`, `dpar` and `nlpar` slots each prior belongs in. It routes on
    to `paired-count-anatomy` for shrinkage priors in the
    interpretable coordinates.
  - *A worked partially observed fit* — the first such example in the package.
    A design with 40% of first counts unrecorded is simulated, fitted, and the
    imputed first count scored against the values withheld from the model.
  - *Migration and errata* — what each release changed, the substitutions for
    existing code, and the two corrected numerical results.

* **The recommended priors in the $(M, f, \delta)$ coordinates have changed.**
  *The anatomy of a paired count* now recommends `normal(0, 0.5)` on `methd`,
  and `normal(0, 1)` with `lb = 0`, a half-normal, on each of `kappas` and
  `kappax`. The previous recommendation was `double_exponential(0, 0.5)` and
  `exponential(1)`. The worked fit in that article was refitted under the new
  priors: every 90% interval still covers its true value, and no posterior
  median moved by more than 1%. `exponential(1)` is documented as an
  alternative of the same shape with a heavier tail.

* **Those priors are no longer called penalised-complexity priors.** They
  follow the same principles: a base model at zero, and a density decreasing
  monotonically away from it. The PC construction additionally requires a
  distance measure derived for the specific model and a constant-rate
  assumption on that distance, neither of which is established for these
  families. The wording is corrected in `README.md`, the get-started vignette,
  *Choosing priors* and *The anatomy of a paired count*.

* **Attributions and citations updated.** The negative-binomial construction is now
  credited to Kirkpatrick & Neale (2016, 2022), and the Poisson case to Holgate
  (1964) and Karlis & Ntzoufras (2003).

* **Precompiled vignette sources renamed.** `<name>.Rmd.orig` is now
  `_<name>.Rmd`, so an editor highlights it as R Markdown; pkgdown skips any
  basename beginning with `_`. The sources do not ship, so this affects
  contributors only.

* **Roxygen markdown is enabled, which changes every help page.** `DESCRIPTION`
  gained `Roxygen: list(markdown = TRUE)`, which this package had never set, so
  roxygen had been passing markdown through to Rd uninterpreted. Backticks
  rendered as literal characters rather than as code formatting,
  cross-references between topics as bracketed text rather than as links, and
  emphasis markers as their own asterisks. Regenerating the six Rd files
  corrected all three.

  One line was restored rather than reformatted. `%` opens a comment in Rd, so
  the unescaped `%in%` in the example that tells the two constructors apart
  truncated its own line at render time:

  ```r
  "vint2" %in% names(brms::standata(fit))   # TRUE for a partially observed fit
  ```

  On `?bipois_partialobs`, `?binegbin_partialobs` and the corresponding
  reference pages of the site, that line ended after `"vint2"`. It is now
  complete.

* **Every documented function states what it returns.** The `bipois` and
  `binegbin` topics document five functions each and had a single `\value{}`
  describing the constructor alone. Both now give the return value and shape of
  `log_lik_*`, `posterior_predict_*` and `posterior_epred_*` separately, as
  well as that of the constructor and its `_stanvars()` companion. The two
  `_partialobs` topics document two functions each, and state that the family
  returned shares its `name` — and therefore its post-processing methods — with
  the fully paired constructor. `Title:` and `Description:` also quote
  `'brms'`, which is the convention CRAN applies to software names.

* **If you hold a fit made with `bipois_cens()`, `binegbin_cens()` or
  `binegbin_joint()`**, install the last version that defined them — 0.9.1 —
  and keep it on the search path for that fit. brms resolves post-processing
  late, by looking up `log_lik_<name>` from the stored family name at call
  time, so `loo()` and `posterior_predict()` need the methods of the version
  the fit was made under. No argument to those calls can route around it, and
  0.10.0 deliberately ships no shim: pinning the whole package is a stronger
  guarantee than a hand-picked subset of forwarders.
# bicountbrms 0.9.1

* **No behaviour change.** Nothing in the four families' likelihoods,
  predictions or expectations differs from 0.9.0. This release adds test
  coverage that was missing and corrects a claim in the README the suite did not
  support. No fitted model is affected and no code needs changing.

* **New tests: prediction on censored rows, with the two excess dispersions
  distinct** (`tests/testthat/test-cens-predict.R`). The `_cens` families exist
  to impute the first margin on rows where it was never observed, and that row
  class had no correctness test. The Monte Carlo check validating
  `posterior_predict` against the exact conditional `P(y1 | y2)` ran at
  `y1_obs == 1` in both `test-binegbin_cens.R` and `test-bipois_cens.R` —
  matched rows, the one class where it cannot fail. Every prediction test also
  built its prep from a single five-dpar `shapex`, which `.SHAPEXONE_NAMES` and
  `.SHAPEXTWO_NAMES` both resolve to, so the two dispersions were never
  distinguished from each other.

  That combination left a specific defect invisible. In
  `posterior_predict_binegbin_cens()`, the conditional split of `N_shared | y2`
  is weighted by `shapextwo` while the fresh `N1` added on top takes
  `shapexone`; exchanging them passed 412 assertions across six test files with
  no failures. The new file fails on it four times. It also pins that `y1_obs`
  does not change the draws — the flag selects a likelihood branch, not a
  prediction — and that `y2 = 0` leaves `posterior_predict` drawing
  `NB2(lambdaone, shapexone)` alone, which is the short-circuit keeping
  `sample()` from being handed an empty support.

* **README: the conditional-prediction identity is now claimed accurately.** The
  Testing section asserted it for the `_cens` families without qualification
  while the suite checked it on matched rows only. It now states that it holds
  on matched *and* censored rows and with the two excess dispersions distinct,
  which is what the tests above pin.

# bicountbrms 0.9.0

* **The package has been split, and renamed.** `pairedcountbrms` carried two
  unrelated suites: families that model a count pair jointly, and families that
  model its difference. They shared no code — not a helper, not a Stan block,
  not a dependency. This release separates them. The joint families are
  `bicountbrms`; the difference families (`skellam1`/`skellam2`,
  `dnorm1`/`dnorm2`, `dlaplace1`/`dlaplace2`) have returned to
  [`skellambrms`](https://github.com/anhsmith/skellambrms), the name under which
  they were originally released.

  **No fitted model needs refitting**, either for the split or for the rename
  below. brms resolves each fit's `log_lik_*` / `posterior_predict_*` /
  `posterior_epred_*` methods off the attached search path at call time, so a
  stored fit works as soon as the package holding its family is attached. No
  dpar name changed. The only source change the split itself requires is
  `library(pairedcountbrms)` → `library(bicountbrms)` or
  `library(skellambrms)`, and any `pairedcountbrms::` prefix.

  The GitHub repository is renamed to `anhsmith/bicountbrms`, and GitHub serves
  a permanent redirect from `anhsmith/pairedcountbrms` for both web and git, so
  existing clones and `pak::pak("anhsmith/pairedcountbrms")` keep working. The
  redirect from `anhsmith/skellambrms` no longer applies: that name is now an
  independent repository holding the difference families, which is what a user
  who installed `skellambrms` was asking for.

  The timing is deliberate. A forthcoming 1.0.0 will be archived and assigned a
  DOI, and cited from a methods article describing the joint families; the
  citable artefact should therefore describe one thing rather than half of two.
  This 0.9.0 makes the split and the additions below available for checking
  first, so that what is archived has been exercised rather than merely tested.

* **Renamed: `binegbin_joint()` → `binegbin_cens()`.** All four families in this
  package model the two counts jointly, so `_joint` named a property common to
  all of them. The property that selects this one is that the first margin may
  be unobserved on some rows, which `_cens` names.
  `binegbin_cens_stanvars()`, `binegbin_cens_lpmf` (Stan) and the three
  post-processing methods follow the same rename.

  **Old code and old fits both keep working.** `binegbin_joint()` and
  `binegbin_joint_stanvars()` remain as deprecated wrappers that warn once per
  call and return the new objects. The three post-processing functions
  `log_lik_binegbin_joint()`, `posterior_predict_binegbin_joint()` and
  `posterior_epred_binegbin_joint()` remain as *silent* forwarders, because a
  fit made before this release stores `family$name == "binegbin_joint"`
  permanently and brms builds its method names from that stored name — nothing
  the caller passes to `loo()` or `posterior_predict()` can route around it. All
  five are removed in the next major version.

  <!-- Amended at 0.10.0: all five were removed there, which is not a major
  version. While the major version is 0, the API is not held out as stable, so
  the removal needs no major bump. Anyone holding such a fit should pin 0.9.1,
  which freezes the whole package rather than a hand-picked subset of
  forwarders. The sentence above is left as written because this file is a
  record of what each release said at the time. -->


  No dpar name changed, so the forwarding is a straight hand-off, and a
  five-dpar pre-0.8.0 fit still resolves through both compatibility layers at
  once (a test pins this).

* **New family: `bipois_cens()` / `bipois_cens_stanvars()`.** The
  censoring-aware bivariate Poisson — `binegbin_cens()`'s equidispersed
  counterpart, and `bipois()`'s censoring-aware one. Three log-linked dpars
  (`mu`, `lambdaone`, `lambdatwo`) and two `vint()` integers,
  `vint(y2, y1_obs)`, exactly as `binegbin_cens()` takes them. It has no
  deprecated aliases: it was added in this release under its final name.

  Its censored branch is closed form. A sum of independent Poissons is Poisson,
  so the `y1`-integrated marginal collapses to `poisson_lpmf(y2 | mu +
  lambdatwo)` — no convolution sum, no cutoff. `binegbin_cens()` cannot do this
  (NB2 + NB2 is not NB2 unless the two components share a success probability),
  which is why it evaluates that branch term by term.

  Two things follow. Users whose counts are equidispersed no longer have to fit
  `binegbin_cens()` with its dispersions pressed against the Poisson boundary,
  where sampling degrades. And the suite gains an independent analytic reference
  for `binegbin_cens()`'s censored branch: driving its three dispersions to
  their Poisson limit must reproduce `bipois_cens()`, checked as a limit — the
  error must shrink in proportion to `1/phi` — rather than at a single
  tolerance. Previously that branch was pinned only by the marginal identity,
  which compares it against the matched branch of the same code and so shares
  any error common to both sums.

* **Breaking (numerical): `posterior_epred_binegbin()` is now exact.** It
  previously substituted the *marginal* shared fraction `mu/(mu + lambdatwo)`
  for the *conditional* one — the `bipois()` answer, which is exact only in the
  Poisson limit — and its docstring described this as a point approximation
  adequate for display. It is not adequate. Over a grid of plausible rates and
  dispersions, with `y2` within one standard deviation of its mean, the
  substitution was in error by more than 5% in about half the settings and more
  than 20% in a third, in either direction depending on which component carried
  the greater dispersion.

  `E[N_shared | y2]` is now evaluated as the mean of the discrete conditional
  over `k = 0..y2` — the same weights `posterior_predict_binegbin()` already
  samples from, summed rather than sampled — so the expectation and the
  predictions agree by construction rather than approximately.

  Numbers produced by an earlier version's `posterior_epred_binegbin()` should
  be recomputed. `log_lik()`, `loo()` and `posterior_predict()` are unaffected;
  the likelihood has not changed.

* **New: `posterior_epred_binegbin_cens()`.** That family previously defined no
  `posterior_epred` at all, on the view that `E[y1]` is ambiguous when `y1` is
  unobserved. It is not: `E[y1 | y2]` is well defined whether or not `y1` was
  recorded, and it is the quantity a user imputing the missing margin wants. It
  is returned for every row, matched and censored alike, matching
  `posterior_predict_binegbin_cens()`'s existing convention. Stored fits reach it
  through the deprecated `posterior_epred_binegbin_joint()` forwarder, so they
  gain the method without refitting.

* **All four families now share one `posterior_epred` convention**, stated once
  in `tests/testthat/test-epred.R` and held to one standard: the expectation
  must equal the mean of that family's own `posterior_predict` draws to within
  Monte Carlo error. Before this release, the three families answered three
  different ways — exactly, approximately, and not at all.

* `skellam` moves from **Imports to Suggests**. No function in `R/` calls it now
  that the Skellam families have left; `tests/testthat/test-bipois.R` still uses
  `skellam::dskellam()` as the independent reference for the induced difference
  distribution, which is worth keeping.

* The coverage gate's environment variable is renamed
  `PAIREDCOUNTBRMS_COVERAGE` → `BICOUNTBRMS_COVERAGE`.

* The README and the getting-started vignette are rewritten around the four
  joint families, and the vignette's "one limitation to know about" section —
  which described brms's `posterior_epred()` failure on *truncated* custom-family
  fits — is replaced by a statement of the epred convention above. That
  limitation cannot arise here: these families do not support `resp_trunc()`. It
  applies to the difference families and is now documented in `skellambrms`.

---

Entries below record this code's history under its former names. It was
`skellambrms` through 0.5.0 and `pairedcountbrms` from 0.6.0 to 0.8.0, and those
packages contained the difference families as well, so some entries describe
families that are no longer here. They are left as written rather than
retrospectively edited. Releases before 0.4.0, which predate the joint families
entirely, are recorded in
[`skellambrms`'s NEWS](https://github.com/anhsmith/skellambrms/blob/master/NEWS.md).

# pairedcountbrms 0.8.0

* **Breaking (dpar split): `binegbin_joint()`'s single excess dispersion
  `shapex` is now the per-margin pair `shapexone`/`shapextwo`**, giving the
  family six dpars (`mu`, `lambdaone`, `lambdatwo`, `shapes`, `shapexone`,
  `shapextwo`), all log-linked. The two source-only excess components may now
  differ in overdispersion, which the single `shapex` forced them not to. The
  rates were already free to differ; the dispersions now are too.

* **The symmetric model is a formula constraint, not a separate family.**
  `shapexone == shapextwo` recovers the pre-0.8.0 likelihood term for term.
  Express it by routing both dpars through one non-linear parameter:

  ```r
  bf(y1 | vint(y2, y1_obs) ~ 1,
     mu ~ 1 + (1 | vessel),
     nlf(shapexone ~ shapexx),
     nlf(shapextwo ~ shapexx),
     shapexx ~ 1, ..., nl = TRUE)
  ```

  This is the only source change a symmetric model needs. Note the prior moves
  with the name: `prior(..., dpar = "shapex")` becomes
  `prior(..., nlpar = "shapexx")`.

* **Fits made before this release still post-process; no refitting and no shim
  are required.** brms resolves a custom family's `log_lik_*` /
  `posterior_predict_*` by name against the live search path at call time, not
  from anything frozen in the fit, so a stored fit always runs the currently
  attached code. Only the dpar NAMES it declares are frozen, and both excess
  dispersions now fall back to a five-dpar fit's single `shapex` — which is
  precisely the constraint that fit was estimated under. Verified on stored
  five-dpar fits: `log_lik()` and `posterior_predict()` agree with 0.7.0
  bitwise, and `loo()` runs unchanged.

* The generalised family absorbs a project-local asymmetric family developed
  outside the package. During the migration, the two were compared element by
  element over 15,480,000 pointwise log-likelihood values from ten stored
  fits, and every value was bitwise identical. That was a one-off check
  against artefacts outside this repository, not a package test.

* `binegbin_mfd_to_dpars()` gains `kappaxone`/`kappaxtwo`, and
  `binegbin_dpars_to_mfd()` gains `shapexone`/`shapextwo`, so the
  $(M, f, \delta)$ converters serve `binegbin_joint()`'s per-margin excess
  dispersions under the family's own dpar names instead of leaving the caller
  to re-derive `shape = 1/kappa^2`. Purely additive: the existing
  `kappax`/`shapex` arguments are unchanged and still return `shapex`/`kappax`.
  Supplying `kappax` together with either of the new arguments is an error,
  since they are two spellings of the same quantity for different families.

* Note on identifiability: `shapextwo` governs the always-observed margin and
  appears on both branches, so the censored rows inform it. `shapexone`
  appears only on the matched branch and is identified solely by the matched
  rows, as `lambdaone` is. A design with few matched rows will learn
  `shapexone` weakly.

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
  `prepare_predictions()` on an older fit returns a `prep` with the previous
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
  credible interval. For a correct model with a calibrated posterior, this is a
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
  `BICOUNTBRMS_COVERAGE=true` (then `PAIREDCOUNTBRMS_COVERAGE`), a named gate rather than a silent default. It
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
  Negative-Binomial model and five dpars as `binegbin()`, but each row supplies
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

