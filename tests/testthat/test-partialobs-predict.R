# tests/testthat/test-partialobs-predict.R
#
# PREDICTION ON UNMATCHED ROWS, AND THE ASYMMETRIC DISPERSIONS.
#
# Two gaps this file closes, both in the prediction path
# rather than in their likelihood.
#
# 1. THE UNMATCHED ROW WAS NEVER CHECKED FOR CORRECTNESS. The Monte Carlo test
#    that validates posterior_predict against the exact conditional
#    P(y1 | y2) = joint / marginal runs at `vint2 = 1L` in both
#    test-binegbin-partialobs.R and test-bipois-partialobs.R. Unmatched rows take the same
#    code path only because posterior_predict_binegbin() ignores y1_obs --
#    true by reading R/binegbin.R, but asserted nowhere, so a later edit
#    that made the function branch on the flag would pass the suite.
#
# 2. EVERY PREDICTION TEST RAN WITH shapexone == shapextwo. test-epred.R and
#    test-binegbin-partialobs.R both build their prep from a single five-dpar
#    `shapex`, which .SHAPEXONE_NAMES and .SHAPEXTWO_NAMES resolve to
#    identically. The two dispersions do different jobs in
#    posterior_predict_binegbin() -- shapextwo weights the conditional
#    split of N_shared | y2, shapexone governs the fresh N1 added on top -- so
#    swapping those two lines passed the entire suite. The lpmf has a dedicated
#    transposition test for exactly this hazard (test-binegbin-dispersions.R,
#    "swapping both rates and both dispersions transposes the joint"); the
#    prediction functions had no equivalent.
#
#    MEASURED, NOT ASSUMED. Injecting that swap into
#    posterior_predict_binegbin() left 412 R-side assertions across
#    test-epred.R, test-binegbin-partialobs.R, test-binegbin-dispersions.R,
#    test-dpar-compat.R, test-deprecated.R and test-bipois-partialobs.R passing, with
#    zero failures. The tests below fail on it four times.
#
# Every check is therefore run on a UNMATCHED row with the two excess
# dispersions an order of magnitude apart, and each is paired with a
# non-vacuity assertion showing it discriminates.
#
# R-side only -- no Stan compilation, so this runs in the fast suite.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2
SX1 <- 0.6    # y1 excess: strongly overdispersed
SX2 <- 7      # y2 excess: mildly overdispersed
Y2  <- 6L
ND  <- 2e5L

asym_prep <- function(vint2, y2 = Y2, sx1 = SX1, sx2 = SX2, ndraws = ND) {
  make_synthetic_prep(
    dpars = list(
      mu        = rep(MU,     ndraws),
      lambdaone = rep(LONE,   ndraws),
      lambdatwo = rep(LTWO,   ndraws),
      shapes    = rep(SHAPES, ndraws),
      shapexone = rep(sx1,    ndraws),
      shapextwo = rep(sx2,    ndraws)
    ),
    Y = 0L, vint1 = y2, vint2 = vint2
  )
}

# Exact P(y1 = x | y2) = joint(x, y2) / marginal(y2), from the R reference.
cond_pmf <- function(xs, y2, sx1, sx2) {
  lp <- binegbin_lpmf_r(xs, rep(y2, length(xs)),
                             MU, LONE, LTWO, SHAPES, sx1, sx2)
  mx <- max(lp)
  exp(lp - (mx + log(sum(exp(lp - mx)))))
}

# ---------------------------------------------------------------------------
# The unmatched row, with the dispersions genuinely different
# ---------------------------------------------------------------------------

test_that("posterior_predict on a UNMATCHED row reproduces the exact conditional", {
  # The row class the existing MC test never exercises, at shapexone far from
  # shapextwo so the two cannot stand in for each other.
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep(vint2 = 0L))
  expect_length(draws, ND)

  K  <- 100
  xs <- 0:K
  p_cond <- cond_pmf(xs, Y2, SX1, SX2)
  emp <- tabulate(draws + 1L, nbins = K + 1L) / ND

  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)
})

test_that("the conditional the draws match is NOT the dispersion-swapped one", {
  # Non-vacuity for the test above, and the check that pins which dispersion
  # enters the conditional weights and which the added N1 draw. Asserted as a
  # property of the two pmfs, so it does not depend on the sampler: if these
  # were close the test above could not discriminate.
  xs <- 0:100
  p_true    <- cond_pmf(xs, Y2, SX1, SX2)
  p_swapped <- cond_pmf(xs, Y2, SX2, SX1)
  expect_gt(max(abs(p_true - p_swapped)), 0.05)
  expect_gt(abs(sum(xs * p_true) - sum(xs * p_swapped)), 0.5)
})

test_that("y1_obs does not change the draws, matched vs unmatched", {
  # posterior_predict imputes y1 on EVERY row; the flag selects a likelihood
  # branch, not a prediction. test-epred.R pins this for the expectation only.
  # Same seed, so this is an identity rather than a distributional claim.
  set.seed(11); a <- posterior_predict_binegbin(1, asym_prep(vint2 = 0L, ndraws = 500L))
  set.seed(11); b <- posterior_predict_binegbin(1, asym_prep(vint2 = 1L, ndraws = 500L))
  expect_identical(a, b)

  pois <- function(v2) make_synthetic_prep(
    dpars = list(mu = rep(MU, 500L), lambdaone = rep(LONE, 500L),
                 lambdatwo = rep(LTWO, 500L)),
    Y = 0L, vint1 = Y2, vint2 = v2)
  set.seed(11); pa <- posterior_predict_bipois(1, pois(0L))
  set.seed(11); pb <- posterior_predict_bipois(1, pois(1L))
  expect_identical(pa, pb)
})

test_that("epred equals the mean of its own draws on a unmatched asymmetric row", {
  # The epred convention (see test-epred.R) applied where it has not been
  # checked: unmatched row, shapexone != shapextwo. posterior_epred reads only
  # shapextwo, so this also pins that it is the SECOND margin's dispersion that
  # enters -- reading .SHAPEXONE_NAMES there would pass every other test.
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep(vint2 = 0L))
  ep    <- unique(as.vector(posterior_epred_binegbin(asym_prep(vint2 = 0L, ndraws = 2L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # Not vacuous: the shapexone-weighted answer is materially different.
  wrong <- unique(as.vector(posterior_epred_binegbin(
    asym_prep(vint2 = 0L, sx1 = SX2, sx2 = SX1, ndraws = 2L))))
  expect_gt(abs(wrong - ep), 0.1)
})

# ---------------------------------------------------------------------------
# y2 = 0, where the shared component is degenerate
# ---------------------------------------------------------------------------

test_that("y2 = 0 makes posterior_predict draw the private component alone", {
  # N_shared | y2 = 0 is degenerate at 0, so the draws must be exactly
  # NB2(lambdaone, shapexone). This is the branch posterior_predict_binegbin()
  # short-circuits (`if (y2 == 0) 0L`), needed because sample(0, 1) would draw
  # from 1:0 -- and it was untested for the prediction function, epred aside.
  #
  # It is also the modal row class in a sparse partially observed design: where y2 = 0 the
  # imputation carries no information from the observed margin at all.
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep(vint2 = 0L, y2 = 0L))
  expect_length(draws, ND)
  expect_true(all(draws >= 0))

  expect_equal(mean(draws), LONE, tolerance = 0.05)
  expect_equal(var(draws), LONE + LONE^2 / SX1, tolerance = 0.5)

  # And the epred agrees exactly, with no Monte Carlo error to absorb.
  ep <- unique(as.vector(posterior_epred_binegbin(
    asym_prep(vint2 = 0L, y2 = 0L, ndraws = 2L))))
  expect_equal(ep, LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# The SAME checks on the one-vint path
# ---------------------------------------------------------------------------
#
# Everything above runs through a two-vint prep. Most users will never build
# one: a fully paired design calls binegbin(), passes vint(y2) alone, and lands
# on a prep with no vint2 at all. That path carries the identical dispersion
# routing and, before 0.10.0, could not be tested at all -- binegbin() had a
# single shapex, so shapexone and shapextwo resolved to the same number and no
# exchange of them was observable.
#
# Two blocks above are deliberately not mirrored here:
#
#   "the conditional ... is NOT the dispersion-swapped one" is a property of
#   the reference pmf alone, not of any method, so it holds for both shapes and
#   is asserted once.
#
#   "y1_obs does not change the draws" has no one-vint counterpart -- there is
#   no flag on this path. Its cross-shape analogue, that the one-vint answer
#   equals the two-vint answer at y1_obs = 1, is test-unified-vint.R's subject.

# The same fixture with no vint2, which is what vint(y2) alone produces.
asym_prep_matched <- function(y2 = Y2, sx1 = SX1, sx2 = SX2, ndraws = ND) {
  asym_prep(vint2 = NULL, y2 = y2, sx1 = sx1, sx2 = sx2, ndraws = ndraws)
}

test_that("posterior_predict on a ONE-VINT row reproduces the exact conditional", {
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep_matched())
  expect_length(draws, ND)

  K  <- 100
  xs <- 0:K
  p_cond <- cond_pmf(xs, Y2, SX1, SX2)
  emp <- tabulate(draws + 1L, nbins = K + 1L) / ND

  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)
})

test_that("epred equals the mean of its own draws on a one-vint asymmetric row", {
  # Pins that posterior_epred_binegbin() reads shapextwo for the conditional
  # split, exactly as its two-vint counterpart does. Before 0.10.0 there was
  # one dispersion here and nothing to get wrong; there are now two.
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep_matched())
  ep    <- unique(as.vector(posterior_epred_binegbin(asym_prep_matched(ndraws = 2L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # Not vacuous: the shapexone-weighted answer is materially different.
  wrong <- unique(as.vector(posterior_epred_binegbin(
    asym_prep_matched(sx1 = SX2, sx2 = SX1, ndraws = 2L))))
  expect_gt(abs(wrong - ep), 0.1)
})

test_that("y2 = 0 on the one-vint path draws the private component alone", {
  set.seed(20260805)
  draws <- posterior_predict_binegbin(1, asym_prep_matched(y2 = 0L))
  expect_length(draws, ND)
  expect_true(all(draws >= 0))

  expect_equal(mean(draws), LONE, tolerance = 0.05)
  expect_equal(var(draws), LONE + LONE^2 / SX1, tolerance = 0.5)

  ep <- unique(as.vector(posterior_epred_binegbin(
    asym_prep_matched(y2 = 0L, ndraws = 2L))))
  expect_equal(ep, LONE, tolerance = 1e-12)
})

test_that("log_lik on the one-vint path routes the two dispersions correctly", {
  # The likelihood counterpart. log_lik_binegbin() must put shapexone on the
  # y1 term and shapextwo on the y2 term; exchanging them changes the value
  # whenever y1 != y2, which is the case here.
  Y1 <- 9L
  prep <- asym_prep_matched(ndraws = 2L)
  prep$data$Y <- Y1

  got <- unique(as.vector(log_lik_binegbin(1, prep)))
  expect_length(got, 1L)

  want <- binegbin_lpmf_r(Y1, Y2, MU, LONE, LTWO, SHAPES, SX1, SX2)
  expect_equal(got, want, tolerance = 1e-10)

  # Not vacuous: the swapped routing gives a materially different value.
  swapped <- binegbin_lpmf_r(Y1, Y2, MU, LONE, LTWO, SHAPES, SX2, SX1)
  expect_gt(abs(swapped - want), 0.05)
})

test_that("bipois predicts and expects the same on one-vint and two-vint preps", {
  # The Poisson families have no dispersion to misroute, so what is at stake
  # here is only that the one-vint prep is handled at all -- the shape nothing
  # tested before 0.10.0.
  pois_prep <- function(v2) make_synthetic_prep(
    dpars = list(mu = rep(MU, 500L), lambdaone = rep(LONE, 500L),
                 lambdatwo = rep(LTWO, 500L)),
    Y = 0L, vint1 = Y2, vint2 = v2)

  set.seed(11); one <- posterior_predict_bipois(1, pois_prep(NULL))
  set.seed(11); two <- posterior_predict_bipois(1, pois_prep(1L))
  expect_identical(one, two)

  ep <- unique(as.vector(posterior_epred_bipois(pois_prep(NULL))))
  expect_equal(ep, Y2 * MU / (MU + LTWO) + LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# NOT YET COVERED -- deliberately recorded rather than left implicit
# ---------------------------------------------------------------------------
#
# 1. THE EPRED IS NOT GRID-CHECKED AGAINST AN INDEPENDENT REFERENCE.
#    posterior_epred_<family>() is deterministic, so it can be compared against
#    a brute-force conditional expectation over a whole parameter grid in
#    seconds -- no Monte Carlo, no sampler. test-epred.R checks it against the
#    mean of posterior_predict draws at ONE setting plus the Poisson limit, and
#    this file adds one unmatched asymmetric setting and y2 = 0. Neither is a
#    grid. The reference to compare against is .e_shared_given_y2_nb()'s
#    definition in R/utils.R, re-implemented independently rather than called.
#    Cheap; the natural next addition to this file.
#
# 2. NO FIT-LEVEL TEST THAT UNMATCHED ROWS DO NOT SHARPEN `shapexone`.
#    test-binegbin-dispersions.R pins at the LIKELIHOOD level that shapexone is
#    absent from the y2-only branch. The predictive counterpart is a property of
#    a fit: as the matched fraction falls, the posteriors for lambdaone and
#    shapexone must WIDEN, and the imputation intervals with them. If adding
#    y2-only rows tightened shapexone, the unmatched branch would be leaking
#    information it cannot have. Needs Stan and several fits, so it belongs
#    behind BICOUNTBRMS_COVERAGE, reusing coverage_recovery()'s compile-once /
#    update(recompile = FALSE) pattern in helper-coverage.R. Still needs a
#    concrete spec: which matched fractions, how many replicates, and what
#    counts as a failure.
#
# 3. PREDICTIVE COVERAGE UNDER CORRECT SPECIFICATION IS DELIBERATELY NOT HERE.
#    Simulating from the generative model, withholding y1, refitting and checking
#    the held-out truth falls inside its predictive interval is nominal BY
#    CONSTRUCTION when the model is correctly specified, so it can only catch
#    implementation bugs that the exact-conditional checks above already catch,
#    far more sharply and about a thousand times faster. Its real value is
#    quantifying how imputation degrades as the matched fraction falls -- which
#    is guidance for a user, not a pass/fail gate, and belongs in a vignette.
#    See the identifiability caveat in README.md.
