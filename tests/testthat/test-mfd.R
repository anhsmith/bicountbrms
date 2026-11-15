# mfd_to_dpars() and dpars_to_mfd(): the map between the congruence-bias
# coordinates (M, f, delta, and kappa) and the dpars of the two families.

test_that("forward map satisfies its defining identities", {
  M <- 12; f <- 0.67; delta <- 0.2
  d <- mfd_to_dpars(M, f, delta)

  # M is the midpoint: mu + (lambdaone + lambdatwo)/2 == M, for any delta
  expect_equal(d$mu + (d$lambdaone + d$lambdatwo) / 2, M)
  # f is the congruence
  expect_equal(d$mu / M, f)
  # delta is the half log ratio of the exclusive rates
  expect_equal(0.5 * log(d$lambdaone / d$lambdatwo), delta)
})

test_that("M stays pinned to the midpoint regardless of bias", {
  # The source bias changes the two exclusive rates but not M.
  for (delta in c(-2, -0.5, 0, 0.5, 2)) {
    d <- mfd_to_dpars(M = 10, f = 0.5, delta = delta)
    expect_equal(d$mu + (d$lambdaone + d$lambdatwo) / 2, 10)
  }
})

test_that("round trip is exact on the interior", {
  grid <- expand.grid(
    M     = c(0.5, 3, 12, 40),
    f     = c(0.05, 0.3, 0.67, 0.95),
    delta = c(-1.5, -0.2, 0, 0.2, 1.5)
  )
  d <- mfd_to_dpars(grid$M, grid$f, grid$delta)
  back <- dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)

  expect_equal(back$M, grid$M)
  expect_equal(back$f, grid$f)
  expect_equal(back$delta, grid$delta)
})

test_that("reverse round trip is exact", {
  mu <- c(8.04, 1, 20); lone <- c(4.75, 0.5, 3); ltwo <- c(3.17, 2, 3)
  m <- dpars_to_mfd(mu, lone, ltwo)
  d <- mfd_to_dpars(m$M, m$f, m$delta)

  expect_equal(d$mu, mu)
  expect_equal(d$lambdaone, lone)
  expect_equal(d$lambdatwo, ltwo)
})

test_that("beta is the fractional imbalance and equals tanh(delta)", {
  m <- dpars_to_mfd(mu = 5, lambdaone = 6, lambdatwo = 2)
  expect_equal(m$beta, (6 - 2) / (6 + 2))
  expect_equal(m$beta, tanh(m$delta))
})

test_that("f = 1 gives zero exclusive rates and an unidentified bias", {
  d <- mfd_to_dpars(M = 12, f = 1, delta = 0.5)
  expect_equal(d$lambdaone, 0)
  expect_equal(d$lambdatwo, 0)
  expect_equal(d$mu, 12)

  # With both exclusive rates at zero, the source bias is undefined, so the
  # inverse returns NA rather than 0.
  back <- dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)
  expect_equal(back$M, 12)
  expect_equal(back$f, 1)
  expect_true(is.na(back$delta))
  expect_true(is.na(back$beta))
})

test_that("f = 0 removes the shared component", {
  d <- mfd_to_dpars(M = 10, f = 0, delta = 0)
  expect_equal(d$mu, 0)
  expect_equal(d$lambdaone, 10)
  expect_equal(d$lambdatwo, 10)
})

test_that("M = 0 leaves f undefined", {
  back <- dpars_to_mfd(mu = 0, lambdaone = 0, lambdatwo = 0)
  expect_equal(back$M, 0)
  expect_true(is.na(back$f))
  expect_true(is.na(back$delta))
})

test_that("one exclusive rate at zero is the finite +/-Inf bias limit", {
  expect_equal(dpars_to_mfd(5, 4, 0)$delta, Inf)
  expect_equal(dpars_to_mfd(5, 0, 4)$delta, -Inf)
  expect_equal(dpars_to_mfd(5, 4, 0)$beta, 1)
  expect_equal(dpars_to_mfd(5, 0, 4)$beta, -1)

  # The forward map with delta = Inf sets lambdatwo to zero.
  d <- mfd_to_dpars(M = 10, f = 0.5, delta = Inf)
  expect_equal(d$lambdatwo, 0)
  expect_equal(d$lambdaone, 10)
})

test_that("dispersion conversion inverts, with the direction reversed", {
  # `kappax` sets both shapexone and shapextwo, because no family has a
  # `shapex` dpar.
  d <- mfd_to_dpars(12, 0.6, 0, kappas = 0.5, kappax = 2)
  expect_equal(d$shapes, 1 / 0.5^2)
  expect_equal(d$shapexone, 1 / 2^2)
  expect_equal(d$shapextwo, 1 / 2^2)
  expect_null(d$shapex)

  back <- dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo,
                       shapes = d$shapes,
                       shapexone = d$shapexone,
                       shapextwo = d$shapextwo)
  expect_equal(back$kappas, 0.5)
  expect_equal(back$kappaxone, 2)
  expect_equal(back$kappaxtwo, 2)

  # raising kappa lowers shape
  expect_lt(mfd_to_dpars(1, 0.5, 0, kappas = 2)$shapes,
            mfd_to_dpars(1, 0.5, 0, kappas = 1)$shapes)
})

test_that("the inverse still reads the single `shapex` in a pre-0.10.0 fit", {
  # dpars_to_mfd() reads the dpars of stored fits, and binegbin fits made
  # before 0.10.0 declare a single `shapex`. The forward direction writes
  # shapexone and shapextwo only.
  back <- dpars_to_mfd(7.2, 2.4, 2.4, shapes = 4, shapex = 1 / 2^2)
  expect_equal(back$kappas, 0.5)
  expect_equal(back$kappax, 2)
  expect_null(back$kappaxone)

  # The same dispersion supplied as shapexone and shapextwo gives the same
  # overdispersion.
  pair <- dpars_to_mfd(7.2, 2.4, 2.4, shapes = 4,
                       shapexone = 1 / 2^2, shapextwo = 1 / 2^2)
  expect_equal(back$kappax, pair$kappaxone)
  expect_equal(back$kappax, pair$kappaxtwo)
})

# --------------------------------------------------------------------------
# Separate dispersions of the two exclusive components
# --------------------------------------------------------------------------

test_that("kappaxone/kappaxtwo convert to shapexone/shapextwo and round-trip", {
  # The converters map kappaxone and kappaxtwo to shapexone and shapextwo by
  # shape = 1/kappa^2, and back.
  d <- mfd_to_dpars(12, 0.6, 0.2, kappas = 0.5,
                    kappaxone = 2, kappaxtwo = 0.25)
  expect_equal(d$shapes,    1 / 0.5^2)
  expect_equal(d$shapexone, 1 / 2^2)
  expect_equal(d$shapextwo, 1 / 0.25^2)
  # No `shapex` is returned.
  expect_null(d$shapex)

  back <- dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo,
                       shapes = d$shapes,
                       shapexone = d$shapexone,
                       shapextwo = d$shapextwo)
  expect_equal(back$kappas,    0.5)
  expect_equal(back$kappaxone, 2)
  expect_equal(back$kappaxtwo, 0.25)
  expect_null(back$kappax)
  # The coordinates M, f, and delta round-trip.
  expect_equal(back$M, 12); expect_equal(back$f, 0.6); expect_equal(back$delta, 0.2)
})

test_that("the pair uses the same kappa <-> shape map as the single kappax", {
  # One value passed as kappax or as kappaxone must give the same shapexone.
  one  <- mfd_to_dpars(12, 0.6, 0, kappax    = 1.7)$shapexone
  pair <- mfd_to_dpars(12, 0.6, 0, kappaxone = 1.7)$shapexone
  expect_identical(one, pair)
  # kappax sets shapexone and shapextwo to the same value.
  tied <- mfd_to_dpars(12, 0.6, 0, kappax = 1.7)
  expect_identical(tied$shapexone, tied$shapextwo)
  expect_equal(mfd_to_dpars(12, 0.6, 0, kappaxone = 0)$shapexone, Inf)
})

test_that("symmetric case: the same kappa twice gives equal dispersions", {
  d <- mfd_to_dpars(12, 0.6, 0, kappaxone = 1.3, kappaxtwo = 1.3)
  expect_identical(d$shapexone, d$shapextwo)
})

test_that("kappax and the per-margin pair cannot be combined", {
  expect_error(
    mfd_to_dpars(12, 0.6, 0, kappax = 1, kappaxone = 2),
    "not both"
  )
  expect_error(
    mfd_to_dpars(12, 0.6, 0, kappax = 1, kappaxtwo = 2),
    "not both"
  )
  expect_error(
    dpars_to_mfd(8, 4, 3, shapex = 1, shapexone = 2),
    "not both"
  )
  # kappax alone and kappaxone alone are accepted.
  expect_type(mfd_to_dpars(12, 0.6, 0, kappax = 1)$shapexone, "double")
  expect_type(mfd_to_dpars(12, 0.6, 0, kappaxone = 1)$shapexone, "double")
})

test_that("the forward direction emits only names a shipping family accepts", {
  # Every name returned by mfd_to_dpars() is a dpar of bipois or binegbin.
  ok <- c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo")
  for (args in list(list(kappas = 0.5, kappax = 1.2),
                    list(kappas = 0.5, kappaxone = 1.2, kappaxtwo = 0.3),
                    list(kappax = 1.2),
                    list())) {
    d <- do.call(mfd_to_dpars, c(list(12, 0.6, 0.1), args))
    expect_true(all(names(d) %in% ok), label = paste(names(d), collapse = ","))
  }
  # The rate names are dpars of bipois(), and all six are dpars of binegbin().
  expect_true(all(setdiff(ok, c("shapes", "shapexone", "shapextwo")) %in% bipois()$dpars))
  expect_true(all(ok %in% binegbin()$dpars))
})

test_that("the forward direction returns dpars in the order the family declares them", {
  # The generated Stan call takes dpars by position, so printed output in a
  # different order would suggest the wrong signature. A chained assignment
  # `a <- b <- value` creates b first, and would place shapextwo before
  # shapexone.
  expect_identical(
    names(mfd_to_dpars(12, 0.67, 0.2, kappas = 0.5, kappax = 1)),
    binegbin()$dpars
  )
  expect_identical(
    names(mfd_to_dpars(12, 0.67, 0.2, kappas = 0.5,
                       kappaxone = 1, kappaxtwo = 0.3)),
    binegbin()$dpars
  )
  expect_identical(
    names(mfd_to_dpars(12, 0.67, 0.2)),
    bipois()$dpars
  )
})

test_that("the pair vectorises over draws like every other argument", {
  d <- mfd_to_dpars(M = c(10, 12), f = 0.6, delta = 0,
                    kappaxone = c(1, 2), kappaxtwo = c(0.5, 0.25))
  expect_length(d$shapexone, 2)
  expect_equal(d$shapexone, 1 / c(1, 2)^2)
  expect_equal(d$shapextwo, 1 / c(0.5, 0.25)^2)
})

test_that("kappa = 0 is the Poisson limit, handled exactly", {
  expect_equal(mfd_to_dpars(12, 0.6, 0, kappas = 0)$shapes, Inf)
  expect_equal(dpars_to_mfd(1, 1, 1, shapes = Inf)$kappas, 0)
})

test_that("arguments recycle, so the map vectorises over draws", {
  d <- mfd_to_dpars(M = c(10, 20, 30), f = 0.5, delta = 0)
  expect_length(d$mu, 3)
  expect_equal(d$mu, c(5, 10, 15))
})

test_that("out-of-range inputs are rejected", {
  expect_error(mfd_to_dpars(M = -1, f = 0.5), "non-negative")
  expect_error(mfd_to_dpars(M = 1, f = 1.5), "\\[0, 1\\]")
  expect_error(dpars_to_mfd(mu = -1, lambdaone = 1, lambdatwo = 1),
               "non-negative")
})

test_that("the map agrees with the trivariate-reduction moment identities", {
  # E[y1] = mu + lambdaone and E[y2] = mu + lambdatwo, so the mean of the two
  # expectations is M, and their difference is lambdaone - lambdatwo.
  d <- mfd_to_dpars(M = 12, f = 0.67, delta = 0.3)
  e_em <- d$mu + d$lambdaone
  e_lb <- d$mu + d$lambdatwo

  expect_equal((e_em + e_lb) / 2, 12)
  expect_equal(e_em - e_lb, d$lambdaone - d$lambdatwo)
})

test_that("the rate half of the map serves bipois unchanged", {
  # bipois has the three rates and no dispersion. The unpaired branch of
  # bipois is Poisson(mu + lambdatwo), so with the converted rates its mean
  # must be E[y2] = mu + lambdatwo = M(1 - (1 - f) tanh(delta)).
  d <- mfd_to_dpars(M = 12, f = 0.67, delta = 0.3)

  ys <- 0:400
  lp <- bipois_lpmf_r(0L, ys, y1_obs = 0L, d$mu, d$lambdaone, d$lambdatwo)
  expect_equal(sum(exp(lp)), 1, tolerance = 1e-10)
  expect_equal(sum(ys * exp(lp)), d$mu + d$lambdatwo, tolerance = 1e-8)
  expect_equal(d$mu + d$lambdatwo, 12 * (1 - (1 - 0.67) * tanh(0.3)))

  # Only the negative-binomial families take a dispersion.
  expect_named(mfd_to_dpars(M = 12, f = 0.67, delta = 0.3),
               c("mu", "lambdaone", "lambdatwo"))
})
