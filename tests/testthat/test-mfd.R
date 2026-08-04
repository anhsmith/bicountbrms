# (M, f, delta) <-> native dpar coordinate transforms.
#
# These are the R-side statement of a map that also appears in the project's
# brms nlf() formulas and in illustrative JS. Testing the round trip and the
# defining identities here is what keeps those in step.

test_that("forward map satisfies its defining identities", {
  M <- 12; f <- 0.67; delta <- 0.2
  d <- binegbin_mfd_to_dpars(M, f, delta)

  # M is the midpoint: mu + (lambdaone + lambdatwo)/2 == M, for any delta
  expect_equal(d$mu + (d$lambdaone + d$lambdatwo) / 2, M)
  # f is the shared share
  expect_equal(d$mu / M, f)
  # delta is the half log ratio of the excesses
  expect_equal(0.5 * log(d$lambdaone / d$lambdatwo), delta)
})

test_that("M stays pinned to the midpoint regardless of bias", {
  # The point of the parameterisation: bias separates the two excess rates but
  # never moves M.
  for (delta in c(-2, -0.5, 0, 0.5, 2)) {
    d <- binegbin_mfd_to_dpars(M = 10, f = 0.5, delta = delta)
    expect_equal(d$mu + (d$lambdaone + d$lambdatwo) / 2, 10)
  }
})

test_that("round trip is exact on the interior", {
  grid <- expand.grid(
    M     = c(0.5, 3, 12, 40),
    f     = c(0.05, 0.3, 0.67, 0.95),
    delta = c(-1.5, -0.2, 0, 0.2, 1.5)
  )
  d <- binegbin_mfd_to_dpars(grid$M, grid$f, grid$delta)
  back <- binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)

  expect_equal(back$M, grid$M)
  expect_equal(back$f, grid$f)
  expect_equal(back$delta, grid$delta)
})

test_that("reverse round trip is exact", {
  mu <- c(8.04, 1, 20); lone <- c(4.75, 0.5, 3); ltwo <- c(3.17, 2, 3)
  m <- binegbin_dpars_to_mfd(mu, lone, ltwo)
  d <- binegbin_mfd_to_dpars(m$M, m$f, m$delta)

  expect_equal(d$mu, mu)
  expect_equal(d$lambdaone, lone)
  expect_equal(d$lambdatwo, ltwo)
})

test_that("beta is the fractional imbalance and equals tanh(delta)", {
  m <- binegbin_dpars_to_mfd(mu = 5, lambdaone = 6, lambdatwo = 2)
  expect_equal(m$beta, (6 - 2) / (6 + 2))
  expect_equal(m$beta, tanh(m$delta))
})

test_that("f = 1 gives zero excesses and an unidentified bias", {
  d <- binegbin_mfd_to_dpars(M = 12, f = 1, delta = 0.5)
  expect_equal(d$lambdaone, 0)
  expect_equal(d$lambdatwo, 0)
  expect_equal(d$mu, 12)

  # Reverse: the bias genuinely cannot be recovered. NA, not 0 -- 0 would
  # assert an unbiased method the data cannot support.
  back <- binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)
  expect_equal(back$M, 12)
  expect_equal(back$f, 1)
  expect_true(is.na(back$delta))
  expect_true(is.na(back$beta))
})

test_that("f = 0 removes the shared component", {
  d <- binegbin_mfd_to_dpars(M = 10, f = 0, delta = 0)
  expect_equal(d$mu, 0)
  expect_equal(d$lambdaone, 10)
  expect_equal(d$lambdatwo, 10)
})

test_that("M = 0 leaves f undefined", {
  back <- binegbin_dpars_to_mfd(mu = 0, lambdaone = 0, lambdatwo = 0)
  expect_equal(back$M, 0)
  expect_true(is.na(back$f))
  expect_true(is.na(back$delta))
})

test_that("one excess at zero is the finite +/-Inf bias limit", {
  expect_equal(binegbin_dpars_to_mfd(5, 4, 0)$delta, Inf)
  expect_equal(binegbin_dpars_to_mfd(5, 0, 4)$delta, -Inf)
  expect_equal(binegbin_dpars_to_mfd(5, 4, 0)$beta, 1)
  expect_equal(binegbin_dpars_to_mfd(5, 0, 4)$beta, -1)

  # ...and the forward map reproduces it from infinite delta
  d <- binegbin_mfd_to_dpars(M = 10, f = 0.5, delta = Inf)
  expect_equal(d$lambdatwo, 0)
  expect_equal(d$lambdaone, 10)
})

test_that("dispersion conversion inverts, with the direction reversed", {
  d <- binegbin_mfd_to_dpars(12, 0.6, 0, kappas = 0.5, kappax = 2)
  expect_equal(d$shapes, 1 / 0.5^2)
  expect_equal(d$shapex, 1 / 2^2)

  back <- binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo,
                                shapes = d$shapes, shapex = d$shapex)
  expect_equal(back$kappas, 0.5)
  expect_equal(back$kappax, 2)

  # raising kappa lowers shape
  expect_lt(binegbin_mfd_to_dpars(1, 0.5, 0, kappas = 2)$shapes,
            binegbin_mfd_to_dpars(1, 0.5, 0, kappas = 1)$shapes)
})

# --------------------------------------------------------------------------
# Per-margin excess dispersions (binegbin_cens, 0.8.0+)
# --------------------------------------------------------------------------

test_that("kappaxone/kappaxtwo convert to shapexone/shapextwo and round-trip", {
  # binegbin_cens carries one excess dispersion per margin, so the converters
  # accept and return the pair under the family's own dpar names rather than
  # making the caller re-derive shape = 1/kappa^2 by hand.
  d <- binegbin_mfd_to_dpars(12, 0.6, 0.2, kappas = 0.5,
                             kappaxone = 2, kappaxtwo = 0.25)
  expect_equal(d$shapes,    1 / 0.5^2)
  expect_equal(d$shapexone, 1 / 2^2)
  expect_equal(d$shapextwo, 1 / 0.25^2)
  # the single-dispersion spelling must NOT appear when the pair was asked for
  expect_null(d$shapex)

  back <- binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo,
                                shapes = d$shapes,
                                shapexone = d$shapexone,
                                shapextwo = d$shapextwo)
  expect_equal(back$kappas,    0.5)
  expect_equal(back$kappaxone, 2)
  expect_equal(back$kappaxtwo, 0.25)
  expect_null(back$kappax)
  # rates round-trip as before
  expect_equal(back$M, 12); expect_equal(back$f, 0.6); expect_equal(back$delta, 0.2)
})

test_that("the pair uses the same kappa <-> shape map as the single kappax", {
  # Not a separate formula: passing one value through kappax and through
  # kappaxone must give the same number, or the two spellings could drift.
  one  <- binegbin_mfd_to_dpars(12, 0.6, 0, kappax    = 1.7)$shapex
  pair <- binegbin_mfd_to_dpars(12, 0.6, 0, kappaxone = 1.7)$shapexone
  expect_identical(one, pair)
  expect_equal(binegbin_mfd_to_dpars(12, 0.6, 0, kappaxone = 0)$shapexone, Inf)
})

test_that("symmetric case: the same kappa twice gives equal dispersions", {
  d <- binegbin_mfd_to_dpars(12, 0.6, 0, kappaxone = 1.3, kappaxtwo = 1.3)
  expect_identical(d$shapexone, d$shapextwo)
})

test_that("kappax and the per-margin pair cannot be combined", {
  expect_error(
    binegbin_mfd_to_dpars(12, 0.6, 0, kappax = 1, kappaxone = 2),
    "not both"
  )
  expect_error(
    binegbin_mfd_to_dpars(12, 0.6, 0, kappax = 1, kappaxtwo = 2),
    "not both"
  )
  expect_error(
    binegbin_dpars_to_mfd(8, 4, 3, shapex = 1, shapexone = 2),
    "not both"
  )
  # ...but either alone is fine
  expect_type(binegbin_mfd_to_dpars(12, 0.6, 0, kappax = 1)$shapex, "double")
  expect_type(binegbin_mfd_to_dpars(12, 0.6, 0, kappaxone = 1)$shapexone, "double")
})

test_that("the pair vectorises over draws like every other argument", {
  d <- binegbin_mfd_to_dpars(M = c(10, 12), f = 0.6, delta = 0,
                             kappaxone = c(1, 2), kappaxtwo = c(0.5, 0.25))
  expect_length(d$shapexone, 2)
  expect_equal(d$shapexone, 1 / c(1, 2)^2)
  expect_equal(d$shapextwo, 1 / c(0.5, 0.25)^2)
})

test_that("kappa = 0 is the Poisson limit, handled exactly", {
  expect_equal(binegbin_mfd_to_dpars(12, 0.6, 0, kappas = 0)$shapes, Inf)
  expect_equal(binegbin_dpars_to_mfd(1, 1, 1, shapes = Inf)$kappas, 0)
})

test_that("arguments recycle, so the map vectorises over draws", {
  d <- binegbin_mfd_to_dpars(M = c(10, 20, 30), f = 0.5, delta = 0)
  expect_length(d$mu, 3)
  expect_equal(d$mu, c(5, 10, 15))
})

test_that("out-of-range inputs are rejected", {
  expect_error(binegbin_mfd_to_dpars(M = -1, f = 0.5), "non-negative")
  expect_error(binegbin_mfd_to_dpars(M = 1, f = 1.5), "\\[0, 1\\]")
  expect_error(binegbin_dpars_to_mfd(mu = -1, lambdaone = 1, lambdatwo = 1),
               "non-negative")
})

test_that("the map agrees with the trivariate-reduction moment identities", {
  # E[y1] = mu + lambdaone and E[y2] = mu + lambdatwo, so the mean of the two
  # expectations is M and their difference is driven entirely by the bias.
  d <- binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.3)
  e_em <- d$mu + d$lambdaone
  e_lb <- d$mu + d$lambdatwo

  expect_equal((e_em + e_lb) / 2, 12)
  expect_equal(e_em - e_lb, d$lambdaone - d$lambdatwo)
})

test_that("the rate half of the map serves bipois_cens unchanged", {
  # bipois_cens takes the same three rates and no dispersion, so only the rate
  # half of the converters applies to it -- there is no kappa to supply. Rather
  # than assert that in prose, feed the converted rates to the family's censored
  # branch and check the resulting distribution has the mean the map predicts:
  # that branch is Poisson(mu + lambdatwo), so its mean must be E[y2] = mu +
  # lambdatwo, which by the identity above is M(1 - (1 - f) * tanh(delta)).
  d <- binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.3)

  ys <- 0:400
  lp <- bipois_cens_lpmf_r(0L, ys, 0L, d$mu, d$lambdaone, d$lambdatwo)
  expect_equal(sum(exp(lp)), 1, tolerance = 1e-10)
  expect_equal(sum(ys * exp(lp)), d$mu + d$lambdatwo, tolerance = 1e-8)
  expect_equal(d$mu + d$lambdatwo, 12 * (1 - (1 - 0.67) * tanh(0.3)))

  # Supplying a dispersion is what distinguishes the negative-binomial
  # families; asking for one here would be a caller error, and there is no
  # bipois-flavoured spelling that would let it through.
  expect_named(binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.3),
               c("mu", "lambdaone", "lambdatwo"))
})
