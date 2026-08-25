# tests/testthat/test-unified-vint.R
#
# THE TWO SHAPES AGREE.
#
# From 0.10.0 each component distribution has one family name and two
# constructors: binegbin() takes vint(y2) and binegbin_partialobs() takes
# vint(y2, y1_obs). Both return the same `name`, so both land on the same lpmf
# and the same three post-processing methods. What distinguishes them at
# post-processing time is nothing but the presence of prep$data$vint2.
#
# The contract this file pins is the one the 0.10.0 work order states as a
# release criterion: for each of log_lik, posterior_predict and
# posterior_epred, the one-vint path on matched data equals the two-vint path
# on the same data with the flag set to 1.
#
# WHY IT NEEDS NON-VACUITY ASSERTIONS. Two of the three methods ignore y1_obs
# by design -- posterior_predict and posterior_epred impute y1 on every row,
# matched or not (see the convention in CLAUDE.md). For those two, "one-vint
# equals two-vint at flag 1" would hold even if vint2 were never read at all,
# so it is asserted alongside the fact that flag 0 gives the same answer too,
# which is the actual specification. log_lik is the method that must branch,
# and there the flag-0 answer must DIFFER -- otherwise the equality above is
# vacuous.
#
# R-side only -- no Stan compilation, so this runs in the fast suite.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2
SX1 <- 0.6    # y1 excess: strongly overdispersed
SX2 <- 7      # y2 excess: mildly overdispersed
Y1  <- 9L     # deliberately != Y2, so a swapped routing is observable
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
# log_lik -- the one method that must read the flag
# ---------------------------------------------------------------------------

test_that("binegbin log_lik: one-vint equals two-vint at y1_obs = 1", {
  one <- unique(as.vector(log_lik_binegbin(1, nb_prep(NULL, ndraws = 2L))))
  two <- unique(as.vector(log_lik_binegbin(1, nb_prep(1L,   ndraws = 2L))))
  expect_length(one, 1L)
  expect_equal(one, two, tolerance = 1e-12)

  # Not vacuous: an absent vint2 must select the MATCHED branch, so the
  # y1_obs = 0 answer has to be a different number. If it were not, the
  # equality above would hold for a method that never read the flag.
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
  # Fixes WHICH branch an absent vint2 selects, against the R reference rather
  # than against the other path. Both families.
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
# posterior_predict -- the flag is ignored, and that IS the specification
# ---------------------------------------------------------------------------

test_that("binegbin posterior_predict: all three prep shapes give one answer", {
  set.seed(101); one  <- posterior_predict_binegbin(1, nb_prep(NULL, ndraws = 500L))
  set.seed(101); two  <- posterior_predict_binegbin(1, nb_prep(1L,   ndraws = 500L))
  set.seed(101); cens <- posterior_predict_binegbin(1, nb_prep(0L,   ndraws = 500L))
  expect_identical(one, two)
  expect_identical(one, cens)
})

test_that("bipois posterior_predict: all three prep shapes give one answer", {
  set.seed(101); one  <- posterior_predict_bipois(1, pois_prep(NULL, ndraws = 500L))
  set.seed(101); two  <- posterior_predict_bipois(1, pois_prep(1L,   ndraws = 500L))
  set.seed(101); cens <- posterior_predict_bipois(1, pois_prep(0L,   ndraws = 500L))
  expect_identical(one, two)
  expect_identical(one, cens)
})

# ---------------------------------------------------------------------------
# posterior_epred -- likewise
# ---------------------------------------------------------------------------

test_that("binegbin posterior_epred: all three prep shapes give one answer", {
  ep <- function(v2) unique(as.vector(posterior_epred_binegbin(nb_prep(v2, ndraws = 2L))))
  expect_equal(ep(NULL), ep(1L), tolerance = 1e-12)
  expect_equal(ep(NULL), ep(0L), tolerance = 1e-12)

  # Not vacuous: the epred is a real conditional expectation, not a constant
  # that any routing would reproduce.
  expect_gt(ep(NULL), LONE)
})

test_that("bipois posterior_epred: all three prep shapes give one answer", {
  ep <- function(v2) unique(as.vector(posterior_epred_bipois(pois_prep(v2, ndraws = 2L))))
  expect_equal(ep(NULL), ep(1L), tolerance = 1e-12)
  expect_equal(ep(NULL), ep(0L), tolerance = 1e-12)
  expect_equal(ep(NULL), Y2 * MU / (MU + LTWO) + LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# The epred convention survives the unification
# ---------------------------------------------------------------------------

test_that("epred equals the mean of posterior_predict draws on the one-vint path", {
  # test-epred.R states this convention for all families; it is repeated here
  # on the prep shape that carries no vint2, with the two dispersions apart.
  set.seed(20260825)
  draws <- posterior_predict_binegbin(1, nb_prep(NULL))
  ep    <- unique(as.vector(posterior_epred_binegbin(nb_prep(NULL, ndraws = 2L))))
  expect_equal(mean(draws), ep, tolerance = 0.1)
})
