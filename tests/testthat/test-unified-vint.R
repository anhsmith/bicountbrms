# tests/testthat/test-unified-vint.R
#
# Agreement between a prep with no vint2, from binegbin() or bipois(), and a
# prep with vint2, from the _partialobs() constructors. The two constructors
# of a family share one name, so brms calls the same methods for both. For
# log_lik, the conditional draw of y1 (.impute_draw_<family>()), and the
# conditional expectation E[y1 | y2] (.impute_mean_<family>()), the prep with
# no vint2 must give the same result as the prep with y1_obs = 1. log_lik
# must give a different result at y1_obs = 0. The conditional draw and
# expectation do not depend on y1_obs, so they must give the same result at
# y1_obs = 0. No Stan compilation is required.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2
SX1 <- 0.6    # y1 exclusive component: strongly overdispersed
SX2 <- 7      # y2 exclusive component: mildly overdispersed
Y1  <- 9L     # differs from Y2, so exchanged dispersions change log_lik
Y2  <- 6L
ND  <- 5e4L

nb_prep <- function(vint2, ndraws = ND, y1 = Y1, y2 = Y2) {
  make_synthetic_prep(
    dpars = list(
      mu        = rep(MU,     ndraws),
      lambdaone = rep(LONE,   ndraws),
      lambdatwo = rep(LTWO,   ndraws),
      shapes    = rep(SHAPES, ndraws),
      shapexone = rep(SX1,    ndraws),
      shapextwo = rep(SX2,    ndraws)
    ),
    Y = y1, vint1 = y2, vint2 = vint2
  )
}

pois_prep <- function(vint2, ndraws = ND, y1 = Y1, y2 = Y2) {
  make_synthetic_prep(
    dpars = list(
      mu        = rep(MU,   ndraws),
      lambdaone = rep(LONE, ndraws),
      lambdatwo = rep(LTWO, ndraws)
    ),
    Y = y1, vint1 = y2, vint2 = vint2
  )
}

# ---------------------------------------------------------------------------
# log_lik, which depends on y1_obs
# ---------------------------------------------------------------------------

test_that("binegbin log_lik: one-vint equals two-vint at y1_obs = 1", {
  one <- unique(as.vector(log_lik_binegbin(1, nb_prep(NULL, ndraws = 2L))))
  two <- unique(as.vector(log_lik_binegbin(1, nb_prep(1L,   ndraws = 2L))))
  expect_length(one, 1L)
  expect_equal(one, two, tolerance = 1e-12)

  # The unpaired branch (y1_obs = 0) differs by more than 0.5, so the
  # equality above would fail for a method that ignored y1_obs.
  cens <- unique(as.vector(log_lik_binegbin(1, nb_prep(0L, ndraws = 2L))))
  expect_gt(abs(cens - one), 0.5)
})

test_that("bipois log_lik: one-vint equals two-vint at y1_obs = 1", {
  one <- unique(as.vector(log_lik_bipois(1, pois_prep(NULL, ndraws = 2L))))
  two <- unique(as.vector(log_lik_bipois(1, pois_prep(1L,   ndraws = 2L))))
  expect_length(one, 1L)
  expect_equal(one, two, tolerance = 1e-12)

  cens <- unique(as.vector(log_lik_bipois(1, pois_prep(0L, ndraws = 2L))))
  expect_gt(abs(cens - one), 0.5)
})

test_that("the one-vint log_lik is the joint, not the marginal", {
  # With no vint2, log_lik must equal the paired branch of the R reference.
  nb <- unique(as.vector(log_lik_binegbin(1, nb_prep(NULL, ndraws = 2L))))
  expect_equal(
    nb,
    binegbin_lpmf_r(Y1, Y2, MU, LONE, LTWO, SHAPES, SX1, SX2),
    tolerance = 1e-10
  )

  po <- unique(as.vector(log_lik_bipois(1, pois_prep(NULL, ndraws = 2L))))
  expect_equal(po, bipois_lpmf_r(Y1, Y2, MU, LONE, LTWO), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Conditional draws, which do not depend on y1_obs
# ---------------------------------------------------------------------------

test_that("binegbin conditional draws: all three prep shapes give one answer", {
  set.seed(101); one  <- .impute_draw_binegbin(1, nb_prep(NULL, ndraws = 500L))
  set.seed(101); two  <- .impute_draw_binegbin(1, nb_prep(1L,   ndraws = 500L))
  set.seed(101); cens <- .impute_draw_binegbin(1, nb_prep(0L,   ndraws = 500L))
  expect_identical(one, two)
  expect_identical(one, cens)
})

test_that("bipois conditional draws: all three prep shapes give one answer", {
  set.seed(101); one  <- .impute_draw_bipois(1, pois_prep(NULL, ndraws = 500L))
  set.seed(101); two  <- .impute_draw_bipois(1, pois_prep(1L,   ndraws = 500L))
  set.seed(101); cens <- .impute_draw_bipois(1, pois_prep(0L,   ndraws = 500L))
  expect_identical(one, two)
  expect_identical(one, cens)
})

# ---------------------------------------------------------------------------
# Conditional expectation, which does not depend on y1_obs
# ---------------------------------------------------------------------------

test_that("binegbin conditional expectation: all three prep shapes give one answer", {
  ep <- function(v2) unique(as.vector(.impute_mean_binegbin(nb_prep(v2, ndraws = 2L))))
  expect_equal(ep(NULL), ep(1L), tolerance = 1e-12)
  expect_equal(ep(NULL), ep(0L), tolerance = 1e-12)

  # The expectation exceeds lambdaone, so it includes the shared component.
  expect_gt(ep(NULL), LONE)
})

test_that("bipois conditional expectation: all three prep shapes give one answer", {
  ep <- function(v2) unique(as.vector(.impute_mean_bipois(pois_prep(v2, ndraws = 2L))))
  expect_equal(ep(NULL), ep(1L), tolerance = 1e-12)
  expect_equal(ep(NULL), ep(0L), tolerance = 1e-12)
  expect_equal(ep(NULL), Y2 * MU / (MU + LTWO) + LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Conditional expectation against the conditional draws
# ---------------------------------------------------------------------------

test_that("the conditional expectation equals the mean of the conditional draws on the one-vint path", {
  # The check in test-impute.R, repeated on a prep with no vint2 and unequal
  # dispersions.
  set.seed(20260825)
  draws <- .impute_draw_binegbin(1, nb_prep(NULL))
  ep    <- unique(as.vector(.impute_mean_binegbin(nb_prep(NULL, ndraws = 2L))))
  expect_equal(mean(draws), ep, tolerance = 0.1)
})
