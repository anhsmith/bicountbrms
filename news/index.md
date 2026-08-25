# Changelog

## bicountbrms 0.9.1

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
  [`posterior_predict_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  the conditional split of `N_shared | y2` is weighted by `shapextwo`
  while the fresh `N1` added on top carries `shapexone`; exchanging them
  passed 412 assertions across six test files with no failures. The new
  file fails on it four times. It also pins that `y1_obs` does not
  change the draws — the flag selects a likelihood branch, not a
  prediction — and that `y2 = 0` leaves `posterior_predict` drawing
  `NB2(lambdaone, shapexone)` alone, which is the short-circuit keeping
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

- **Renamed:
  [`binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  →
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md).**
  All four families in this package model the two counts jointly, so
  `_joint` named a property common to all of them. The property that
  selects this one is that the first margin may be unobserved on some
  rows, which `_cens` names.
  [`binegbin_cens_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md),
  `binegbin_cens_lpmf` (Stan) and the three post-processing methods
  follow the same rename.

  **Old code and old fits both keep working.**
  [`binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  and
  [`binegbin_joint_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  remain as deprecated wrappers that warn once per call and return the
  new objects. The three post-processing functions
  [`log_lik_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md),
  [`posterior_predict_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  and
  [`posterior_epred_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
  remain as *silent* forwarders, because a fit made before this release
  stores `family$name == "binegbin_joint"` permanently and brms builds
  its method names from that stored name — nothing the caller passes to
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) or
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  can route around it. All five are removed in the next major version.

  No dpar name changed, so the forwarding is a straight hand-off, and a
  five-dpar pre-0.8.0 fit still resolves through both compatibility
  layers at once (a test pins this).

- **New family:
  [`bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md)
  /
  [`bipois_cens_stanvars()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md).**
  The censoring-aware bivariate Poisson —
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)’s
  equidispersed counterpart, and
  [`bipois()`](https://anhsmith.github.io/bicountbrms/reference/bipois.md)’s
  censoring-aware one. Three log-linked dpars (`mu`, `lambdaone`,
  `lambdatwo`) and two `vint()` integers, `vint(y2, y1_obs)`, exactly as
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  takes them. It has no deprecated aliases: it was added in this release
  under its final name.

  Its censored branch is closed form. A sum of independent Poissons is
  Poisson, so the `y1`-integrated marginal collapses to
  `poisson_lpmf(y2 | mu + lambdatwo)` — no convolution sum, no cutoff.
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  cannot do this (NB2 + NB2 is not NB2 unless the two components share a
  success probability), which is why it evaluates that branch term by
  term.

  Two things follow. Users whose counts are equidispersed no longer have
  to fit
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)
  with its dispersions pressed against the Poisson boundary, where
  sampling degrades. And the suite gains an independent analytic
  reference for
  [`binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)’s
  censored branch: driving its three dispersions to their Poisson limit
  must reproduce
  [`bipois_cens()`](https://anhsmith.github.io/bicountbrms/reference/bipois_cens.md),
  checked as a limit — the error must shrink in proportion to `1/phi` —
  rather than at a single tolerance. Previously that branch was pinned
  only by the marginal identity, which compares it against the matched
  branch of the same code and so shares any error common to both sums.

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

- **New:
  [`posterior_epred_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md).**
  That family previously defined no `posterior_epred` at all, on the
  view that `E[y1]` is ambiguous when `y1` is unobserved. It is not:
  `E[y1 | y2]` is well defined whether or not `y1` was recorded, and it
  is the quantity a user imputing the missing margin wants. It is
  returned for every row, matched and censored alike, matching
  [`posterior_predict_binegbin_cens()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_cens.md)’s
  existing convention. Stored fits reach it through the deprecated
  [`posterior_epred_binegbin_joint()`](https://anhsmith.github.io/bicountbrms/reference/binegbin_joint-deprecated.md)
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
  It applies to the difference families and now lives in `skellambrms`.

------------------------------------------------------------------------

Entries below record this code’s history under its former names. It was
`skellambrms` through 0.5.0 and `pairedcountbrms` from 0.6.0 to 0.8.0,
and those packages contained the difference families as well, so some
entries describe families that are no longer here. They are left as
written rather than retrospectively edited. Releases before 0.4.0, which
predate the joint families entirely, are recorded in [`skellambrms`’s
NEWS](https://github.com/anhsmith/skellambrms/blob/master/NEWS.md).
