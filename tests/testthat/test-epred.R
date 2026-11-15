# tests/testthat/test-epred.R
#
# posterior_epred_<family>() returns the marginal mean of the first count,
#
#   E[y1] = mu + lambdaone,
#
# on every row, and posterior_predict_<family>() draws y1 from the marginal
# distribution of y1. Neither method reads the recorded y2 or the observation
# indicator. The conditional mean E[y1 | y2] is tested in test-impute.R.
#
# For each family, the epred must equal the mean of a large posterior_predict
# sample, and the variance of the sample must equal Var(N_shared) + Var(N_1).
# The source-exclusive dispersions differ by an order of magnitude, so a draw
# that read shapextwo in place of shapexone would fail. No Stan compilation
# is required.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2
SX1 <- 0.6    # y1 exclusive component: strongly overdispersed
SX2 <- 7      # y2 exclusive component: mildly overdispersed
Y2 <- 9L
ND <- 1e5L    # MC standard error on the mean is ~0.02 at these settings

pois_prep <- function(vint2 = NULL, y2 = Y2, ndraws = ND) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ndraws), lambdaone = rep(LONE, ndraws),
                 lambdatwo = rep(LTWO, ndraws)),
    Y = 0L, vint1 = y2, vint2 = vint2
  )
}

nb_prep <- function(vint2 = NULL, y2 = Y2, ndraws = ND) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ndraws), lambdaone = rep(LONE, ndraws),
                 lambdatwo = rep(LTWO, ndraws), shapes = rep(SHAPES, ndraws),
                 shapexone = rep(SX1, ndraws), shapextwo = rep(SX2, ndraws)),
    Y = 0L, vint1 = y2, vint2 = vint2
  )
}

# ---------------------------------------------------------------------------
# Both families have one
# ---------------------------------------------------------------------------

test_that("both joint families export a posterior_epred method", {
  # brms finds posterior_epred_<name>() by name, and a missing method errors
  # only when called. Each family has one name, shared by its two
  # constructors.
  for (fam in c("bipois", "binegbin")) {
    expect_true(
      is.function(get0(paste0("posterior_epred_", fam),
                       envir = asNamespace("bicountbrms"))),
      label = paste0("posterior_epred_", fam)
    )
  }

  # The _partialobs() constructors return the family names "bipois" and
  # "binegbin", and no methods exist under the removed names.
  expect_identical(bipois_partialobs()$name,   "bipois")
  expect_identical(binegbin_partialobs()$name, "binegbin")
  for (fam in c("bipois_cens", "binegbin_cens", "binegbin_joint")) {
    expect_null(get0(paste0("posterior_epred_", fam),
                     envir = asNamespace("bicountbrms")),
                label = paste0("posterior_epred_", fam, " is gone"))
  }
})

# ---------------------------------------------------------------------------
# The marginal mean, exactly, under every prep shape
# ---------------------------------------------------------------------------

test_that("the epred of both families is mu + lambdaone under every prep shape", {
  for (v2 in list(NULL, 1L, 0L)) {
    expect_equal(unique(as.vector(posterior_epred_bipois(pois_prep(v2, ndraws = 2L)))),
                 MU + LONE, tolerance = 1e-12)
    expect_equal(unique(as.vector(posterior_epred_binegbin(nb_prep(v2, ndraws = 2L)))),
                 MU + LONE, tolerance = 1e-12)
  }
})

test_that("the epred does not depend on the recorded y2", {
  # The conditional expectation changes by more than 1 between y2 = 0 and
  # y2 = 20, so a method that conditioned on y2 would fail.
  for (y2 in c(0L, 3L, 20L)) {
    expect_equal(posterior_epred_bipois(pois_prep(y2 = y2, ndraws = 2L)),
                 posterior_epred_bipois(pois_prep(ndraws = 2L)))
    expect_equal(posterior_epred_binegbin(nb_prep(y2 = y2, ndraws = 2L)),
                 posterior_epred_binegbin(nb_prep(ndraws = 2L)))
  }
  expect_gt(
    abs(.impute_mean_binegbin(nb_prep(y2 = 20L, ndraws = 2L))[1, 1] -
          .impute_mean_binegbin(nb_prep(y2 = 0L, ndraws = 2L))[1, 1]),
    1
  )
})

test_that("posterior_predict does not read the recorded y2 or the indicator", {
  # The seed is the same for both calls, so the draws must be identical.
  set.seed(11); a <- posterior_predict_binegbin(1, nb_prep(NULL, y2 = 0L,  ndraws = 500L))
  set.seed(11); b <- posterior_predict_binegbin(1, nb_prep(0L,   y2 = 20L, ndraws = 500L))
  expect_identical(a, b)
  set.seed(11); a <- posterior_predict_bipois(1, pois_prep(NULL, y2 = 0L,  ndraws = 500L))
  set.seed(11); b <- posterior_predict_bipois(1, pois_prep(0L,   y2 = 20L, ndraws = 500L))
  expect_identical(a, b)
})

# ---------------------------------------------------------------------------
# epred equals the mean of the draws, and the draws have the marginal variance
# ---------------------------------------------------------------------------

test_that("the bipois epred and variance match its posterior_predict draws", {
  set.seed(20261114)
  draws <- posterior_predict_bipois(1, pois_prep())
  ep    <- unique(as.vector(posterior_epred_bipois(pois_prep(ndraws = 2L))))
  expect_equal(mean(draws), ep, tolerance = 0.01)
  # y1 ~ Poisson(mu + lambdaone), so Var(y1) = E[y1].
  expect_equal(var(draws), MU + LONE, tolerance = 0.03)
})

test_that("the binegbin epred and variance match its posterior_predict draws", {
  set.seed(20261114)
  draws <- posterior_predict_binegbin(1, nb_prep())
  ep    <- unique(as.vector(posterior_epred_binegbin(nb_prep(ndraws = 2L))))
  expect_equal(mean(draws), ep, tolerance = 0.01)

  # Var(y1) = Var(N_shared) + Var(N1), each NB2 variance m + m^2 / phi.
  var_y1 <- MU + MU^2 / SHAPES + LONE + LONE^2 / SX1
  expect_equal(var(draws), var_y1, tolerance = 0.05)

  # Reading shapextwo for N_1 would give a variance of 21.8 rather than 35.5,
  # far outside the tolerance.
  var_swapped <- MU + MU^2 / SHAPES + LONE + LONE^2 / SX2
  expect_gt(abs(var_swapped - var_y1) / var_y1, 0.3)
})

test_that("the binegbin marginal draws reduce to the Poisson ones as the dispersions grow", {
  # NB2(m, phi) -> Poisson(m) as phi -> Inf, so the marginal variance of the
  # binegbin draw must approach mu + lambdaone.
  p <- make_synthetic_prep(
    dpars = list(mu = rep(MU, ND), lambdaone = rep(LONE, ND),
                 lambdatwo = rep(LTWO, ND), shapes = rep(1e7, ND),
                 shapexone = rep(1e7, ND), shapextwo = rep(1e7, ND)),
    Y = 0L, vint1 = Y2
  )
  set.seed(20261114)
  draws <- posterior_predict_binegbin(1, p)
  expect_equal(var(draws), MU + LONE, tolerance = 0.03)
})
