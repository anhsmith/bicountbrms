# tests/testthat/test-partialobs-predict.R
#
# Conditional prediction of y1 given y2 on unpaired rows, with unequal
# dispersions of the two source-exclusive components. In
# .impute_draw_binegbin(), shapextwo enters the conditional distribution of
# N_shared given y2, and shapexone enters the new draw of N_1. The tests use
# dispersions an order of magnitude apart, so an exchange of the two would
# fail, and each test has a check that it is not vacuous. No Stan compilation
# is required.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2
SX1 <- 0.6    # y1 exclusive component: strongly overdispersed
SX2 <- 7      # y2 exclusive component: mildly overdispersed
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
# Unpaired rows with unequal dispersions
# ---------------------------------------------------------------------------

test_that("the conditional draw on an UNPAIRED row reproduces the exact conditional", {
  # The empirical distribution of 2e5 draws on an unpaired row must match the
  # exact conditional pmf.
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep(vint2 = 0L))
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
  # The conditional pmf with the two dispersions exchanged differs from the
  # correct pmf by more than the tolerance of the test above, so that test
  # detects an exchange.
  xs <- 0:100
  p_true    <- cond_pmf(xs, Y2, SX1, SX2)
  p_swapped <- cond_pmf(xs, Y2, SX2, SX1)
  expect_gt(max(abs(p_true - p_swapped)), 0.05)
  expect_gt(abs(sum(xs * p_true) - sum(xs * p_swapped)), 0.5)
})

test_that("y1_obs does not change the draws, paired vs unpaired", {
  # The observation indicator selects a branch of the likelihood and does not
  # alter the conditional draw. The seed is the same for both calls, so the
  # draws must be identical. test-unified-vint.R checks the expectation.
  set.seed(11); a <- .impute_draw_binegbin(1, asym_prep(vint2 = 0L, ndraws = 500L))
  set.seed(11); b <- .impute_draw_binegbin(1, asym_prep(vint2 = 1L, ndraws = 500L))
  expect_identical(a, b)

  pois <- function(v2) make_synthetic_prep(
    dpars = list(mu = rep(MU, 500L), lambdaone = rep(LONE, 500L),
                 lambdatwo = rep(LTWO, 500L)),
    Y = 0L, vint1 = Y2, vint2 = v2)
  set.seed(11); pa <- .impute_draw_bipois(1, pois(0L))
  set.seed(11); pb <- .impute_draw_bipois(1, pois(1L))
  expect_identical(pa, pb)
})

test_that("epred equals the mean of its own draws on an unpaired asymmetric row", {
  # The conditional expectation must equal the mean of the conditional draws
  # on an unpaired row with unequal dispersions. .impute_mean_binegbin()
  # reads shapextwo alone.
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep(vint2 = 0L))
  ep    <- unique(as.vector(.impute_mean_binegbin(asym_prep(vint2 = 0L, ndraws = 2L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # With the dispersions exchanged, the expectation differs by more than 0.1.
  wrong <- unique(as.vector(.impute_mean_binegbin(
    asym_prep(vint2 = 0L, sx1 = SX2, sx2 = SX1, ndraws = 2L))))
  expect_gt(abs(wrong - ep), 0.1)
})

# ---------------------------------------------------------------------------
# y2 = 0, where the shared component is degenerate
# ---------------------------------------------------------------------------

test_that("y2 = 0 reduces the conditional draw to the exclusive component", {
  # N_shared given y2 = 0 is 0, so the draws follow NB2(lambdaone, shapexone).
  # .impute_draw_binegbin() treats y2 == 0 separately, because sample(0, 1)
  # would draw from 1:0. The tolerances allow for Monte Carlo error.
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep(vint2 = 0L, y2 = 0L))
  expect_length(draws, ND)
  expect_true(all(draws >= 0))

  expect_equal(mean(draws), LONE, tolerance = 0.05)
  expect_equal(var(draws), LONE + LONE^2 / SX1, tolerance = 0.5)

  # The conditional expectation is lambdaone exactly.
  ep <- unique(as.vector(.impute_mean_binegbin(
    asym_prep(vint2 = 0L, y2 = 0L, ndraws = 2L))))
  expect_equal(ep, LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# The same checks on a prep with no vint2
# ---------------------------------------------------------------------------
#
# A fit from binegbin() passes vint(y2) alone, so its prep has no vint2. The
# check on the exchanged pmf depends on the reference alone and is not
# repeated. A prep with no vint2 has no observation indicator; agreement with
# a prep at y1_obs = 1 is tested in test-unified-vint.R.

# The same fixture with no vint2, as vint(y2) alone produces.
asym_prep_paired <- function(y2 = Y2, sx1 = SX1, sx2 = SX2, ndraws = ND) {
  asym_prep(vint2 = NULL, y2 = y2, sx1 = sx1, sx2 = sx2, ndraws = ndraws)
}

test_that("the conditional draw on a ONE-VINT row reproduces the exact conditional", {
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep_paired())
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
  # .impute_mean_binegbin() must use shapextwo for the split of y2 into shared
  # and exclusive counts when the prep has no vint2.
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep_paired())
  ep    <- unique(as.vector(.impute_mean_binegbin(asym_prep_paired(ndraws = 2L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # With the dispersions exchanged, the expectation differs by more than 0.1.
  wrong <- unique(as.vector(.impute_mean_binegbin(
    asym_prep_paired(sx1 = SX2, sx2 = SX1, ndraws = 2L))))
  expect_gt(abs(wrong - ep), 0.1)
})

test_that("y2 = 0 on the one-vint path draws only the exclusive component", {
  set.seed(20260805)
  draws <- .impute_draw_binegbin(1, asym_prep_paired(y2 = 0L))
  expect_length(draws, ND)
  expect_true(all(draws >= 0))

  expect_equal(mean(draws), LONE, tolerance = 0.05)
  expect_equal(var(draws), LONE + LONE^2 / SX1, tolerance = 0.5)

  ep <- unique(as.vector(.impute_mean_binegbin(
    asym_prep_paired(y2 = 0L, ndraws = 2L))))
  expect_equal(ep, LONE, tolerance = 1e-12)
})

test_that("log_lik on the one-vint path routes the two dispersions correctly", {
  # log_lik_binegbin() must apply shapexone to the y1 term and shapextwo to
  # the y2 term. With y1 = 9 and y2 = 6, an exchange changes the value.
  Y1 <- 9L
  prep <- asym_prep_paired(ndraws = 2L)
  prep$data$Y <- Y1

  got <- unique(as.vector(log_lik_binegbin(1, prep)))
  expect_length(got, 1L)

  want <- binegbin_lpmf_r(Y1, Y2, MU, LONE, LTWO, SHAPES, SX1, SX2)
  expect_equal(got, want, tolerance = 1e-10)

  # The exchanged dispersions change the value by more than 0.05.
  swapped <- binegbin_lpmf_r(Y1, Y2, MU, LONE, LTWO, SHAPES, SX2, SX1)
  expect_gt(abs(swapped - want), 0.05)
})

test_that("bipois predicts and expects the same on one-vint and two-vint preps", {
  # bipois has no dispersion, so the test checks only that a prep with no
  # vint2 gives the same draws and expectation as a paired prep.
  pois_prep <- function(v2) make_synthetic_prep(
    dpars = list(mu = rep(MU, 500L), lambdaone = rep(LONE, 500L),
                 lambdatwo = rep(LTWO, 500L)),
    Y = 0L, vint1 = Y2, vint2 = v2)

  set.seed(11); one <- .impute_draw_bipois(1, pois_prep(NULL))
  set.seed(11); two <- .impute_draw_bipois(1, pois_prep(1L))
  expect_identical(one, two)

  ep <- unique(as.vector(.impute_mean_bipois(pois_prep(NULL))))
  expect_equal(ep, Y2 * MU / (MU + LTWO) + LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Not tested
# ---------------------------------------------------------------------------
#
# (i) .impute_mean_binegbin() is not compared with a brute-force conditional
# expectation over a parameter grid. (ii) No fitted-model test checks that
# the posterior of shapexone widens as the proportion of paired rows falls.
