# tests/testthat/test-epred.R
#
# THE EPRED CONVENTION, STATED ONCE FOR ALL FOUR FAMILIES.
#
# Every joint family returns the same quantity from posterior_epred_<family>():
#
#   E[y1 | y2] = E[N_shared | y2] + lambdaone
#
# -- the conditional expectation of the first margin given the observed second
# one, exactly, for every row including the censored ones. Before 0.9.0 the
# three families answered three different ways: bipois exactly, binegbin with
# the marginal shared fraction substituted for the conditional one, and
# binegbin_cens not at all. This file pins the convention so they cannot drift
# apart again.
#
# The standard applied to each is the same and does not depend on knowing the
# closed form: the epred must equal the mean of a large posterior_predict
# sample, because posterior_predict simulates precisely the distribution the
# epred is the mean of. Where a closed form does exist (the Poisson families),
# it is checked against that too, to machine precision rather than to Monte
# Carlo error.
#
# R-side only -- no Stan compilation, so this runs in the fast suite.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2; SHAPEX <- 1.5
Y2 <- 9L
ND <- 5e4L   # MC standard error on the mean is ~0.02 at these settings

pois_prep <- function(vint2 = NULL) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ND), lambdaone = rep(LONE, ND),
                 lambdatwo = rep(LTWO, ND)),
    Y = 0L, vint1 = Y2, vint2 = vint2
  )
}

nb_prep <- function(vint2 = NULL, shapex = SHAPEX) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ND), lambdaone = rep(LONE, ND),
                 lambdatwo = rep(LTWO, ND), shapes = rep(SHAPES, ND),
                 shapex = rep(shapex, ND)),
    Y = 0L, vint1 = Y2, vint2 = vint2
  )
}

# ---------------------------------------------------------------------------
# Every family has one
# ---------------------------------------------------------------------------

test_that("all four joint families export a posterior_epred method", {
  # brms locates these by name convention, so their existence IS the interface.
  # A family without one silently falls back to brms's own dispatch and errors
  # at call time rather than at load time, which is why this is asserted
  # directly rather than left to the end-to-end tests.
  for (fam in c("bipois", "bipois_cens", "binegbin", "binegbin_cens")) {
    expect_true(
      is.function(get0(paste0("posterior_epred_", fam),
                       envir = asNamespace("bicountbrms"))),
      label = paste0("posterior_epred_", fam)
    )
  }
})

# ---------------------------------------------------------------------------
# The Poisson families: exact, and equal to the closed form
# ---------------------------------------------------------------------------

test_that("bipois and bipois_cens epreds equal the Binomial-split closed form", {
  # Conditioning a sum of independent Poissons on its total gives a Binomial,
  # so E[N_shared | y2] = y2 * mu/(mu + lambdatwo) with no sum to evaluate.
  closed <- Y2 * MU / (MU + LTWO) + LONE

  expect_equal(unique(as.vector(posterior_epred_bipois(pois_prep()))),
               closed, tolerance = 1e-12)
  expect_equal(unique(as.vector(posterior_epred_bipois_cens(pois_prep(1L)))),
               closed, tolerance = 1e-12)

  # And the two agree with each other, matched branch or censored: bipois_cens
  # is bipois with the censoring flag attached, and the conditional expectation
  # does not depend on that flag.
  expect_equal(posterior_epred_bipois_cens(pois_prep(1L)),
               posterior_epred_bipois_cens(pois_prep(0L)))
  expect_equal(posterior_epred_bipois(pois_prep()),
               posterior_epred_bipois_cens(pois_prep(1L)))
})

# ---------------------------------------------------------------------------
# The negative-binomial families: exact, checked against their own predictions
# ---------------------------------------------------------------------------

test_that("binegbin epred equals the mean of its posterior_predict draws", {
  # The check that does not presume a closed form. Until 0.9.0 this failed:
  # the epred used the MARGINAL shared fraction mu/(mu + lambdatwo) where the
  # CONDITIONAL one belongs, which at these settings is off by ~0.6% and, over
  # a grid of plausible rates and dispersions, by more than 5% in about half of
  # them and more than 20% in a third -- in either direction, depending on
  # which component carries the greater dispersion.
  set.seed(20260804)
  draws <- posterior_predict_binegbin(1, nb_prep())
  ep    <- unique(as.vector(posterior_epred_binegbin(nb_prep())))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # The superseded approximation, for the record: it must NOT reproduce the
  # exact value, or this test would be vacuous.
  approx <- Y2 * MU / (MU + LTWO) + LONE
  expect_gt(abs(approx - ep), 0.01)
})

test_that("binegbin_cens epred equals the mean of its posterior_predict draws", {
  set.seed(20260804)
  draws <- posterior_predict_binegbin_cens(1, nb_prep(vint2 = 1L))
  ep    <- unique(as.vector(posterior_epred_binegbin_cens(nb_prep(vint2 = 1L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)
})

test_that("the censoring flag does not change either joint family's epred", {
  # posterior_predict imputes y1 on every row, censored included; epred must
  # match that convention, or the two become incomparable row by row exactly
  # where imputation matters.
  expect_equal(posterior_epred_binegbin_cens(nb_prep(vint2 = 1L)),
               posterior_epred_binegbin_cens(nb_prep(vint2 = 0L)))
  expect_equal(posterior_epred_bipois_cens(pois_prep(1L)),
               posterior_epred_bipois_cens(pois_prep(0L)))
})

test_that("binegbin epred equals binegbin_cens epred on a matched row", {
  # The families share a matched-branch likelihood, so they must share a
  # conditional expectation. binegbin_cens resolves its second-margin
  # dispersion through .SHAPEXTWO_NAMES, whose `shapex` fallback is what makes
  # the five-dpar prep here serve both.
  expect_equal(posterior_epred_binegbin(nb_prep()),
               posterior_epred_binegbin_cens(nb_prep(vint2 = 1L)))
})

# ---------------------------------------------------------------------------
# The families agree in the limit where they should
# ---------------------------------------------------------------------------

test_that("the negative-binomial epreds reduce to the Poisson ones as shapes grow", {
  # NB2(m, phi) -> Poisson(m) as phi -> Inf, so binegbin's exact conditional
  # expectation must approach bipois's Binomial closed form. This would NOT
  # have held for the superseded approximation, which equalled the bipois
  # answer at every phi -- it was the Poisson-limit value used everywhere.
  ref <- posterior_epred_bipois(pois_prep())

  err <- function(phi) {
    p <- make_synthetic_prep(
      dpars = list(mu = rep(MU, 2L), lambdaone = rep(LONE, 2L),
                   lambdatwo = rep(LTWO, 2L), shapes = rep(phi, 2L),
                   shapex = rep(phi, 2L)),
      Y = 0L, vint1 = Y2
    )
    max(abs(posterior_epred_binegbin(p) - ref[1:2, , drop = FALSE]))
  }

  e5 <- err(1e5)
  e7 <- err(1e7)
  expect_lt(e5, 1e-2)
  expect_lt(e7, 1e-4)
  expect_lt(e7, e5)          # a limit, not a coincidence
})

test_that("y2 = 0 leaves the conditional expectation at lambdaone alone", {
  # N_shared | y2 = 0 is degenerate at 0 for every family, so the epred reduces
  # to lambdaone exactly. This is the branch .e_shared_given_y2_nb() short-
  # circuits, and getting it wrong would be invisible in the grid checks above.
  zero_pois <- make_synthetic_prep(
    dpars = list(mu = rep(MU, 2L), lambdaone = rep(LONE, 2L),
                 lambdatwo = rep(LTWO, 2L)),
    Y = 0L, vint1 = 0L, vint2 = 1L
  )
  zero_nb <- make_synthetic_prep(
    dpars = list(mu = rep(MU, 2L), lambdaone = rep(LONE, 2L),
                 lambdatwo = rep(LTWO, 2L), shapes = rep(SHAPES, 2L),
                 shapex = rep(SHAPEX, 2L)),
    Y = 0L, vint1 = 0L, vint2 = 1L
  )
  for (ep in list(posterior_epred_bipois(zero_pois),
                  posterior_epred_bipois_cens(zero_pois),
                  posterior_epred_binegbin(zero_nb),
                  posterior_epred_binegbin_cens(zero_nb))) {
    expect_equal(unique(as.vector(ep)), LONE, tolerance = 1e-12)
  }
})
