# Resolution of old dpar names in stored fits. Fits made before 0.7.0 name the
# exclusive rates lambdaem and lambdalb rather than lambdaone and lambdatwo,
# and .get_rate() resolves either spelling. These tests check that the two
# spellings give identical output, and that a prep with neither spelling stops
# with an error naming the dpars found. They use make_synthetic_prep()
# (helper-prep.R) and run without Stan.

# Length-3 vectors give 3 posterior draws, on 4 observations.
v <- list(mu = c(6, 6.5, 5.5), one = c(3, 3.2, 2.8), two = c(2, 2.1, 1.9),
          shapes = c(4, 4.5, 3.5), shapex = c(5, 5.5, 4.5))
Y  <- c(5L, 9L, 0L, 14L)
V1 <- c(4L, 7L, 2L, 11L)
V2 <- c(1L, 0L, 1L, 0L)

# The same values under either spelling; `nm` gives the two rate names.
mk <- function(nm, five = TRUE, vint2 = NULL) {
  dp <- list(mu = v$mu)
  dp[[nm[1]]] <- v$one
  dp[[nm[2]]] <- v$two
  if (five) {
    dp$shapes <- v$shapes
    dp$shapex <- v$shapex
  }
  make_synthetic_prep(dpars = dp, Y = Y, vint1 = V1, vint2 = vint2)
}

NEW <- c("lambdaone", "lambdatwo")
OLD <- c("lambdaem", "lambdalb")

test_that("log_lik is identical under old and new rate-dpar names", {
  for (i in seq_along(Y)) {
    expect_equal(
      log_lik_binegbin(i, mk(OLD)),
      log_lik_binegbin(i, mk(NEW)),
      label = paste0("binegbin, obs ", i)
    )
    expect_equal(
      log_lik_bipois(i, mk(OLD, five = FALSE)),
      log_lik_bipois(i, mk(NEW, five = FALSE)),
      label = paste0("bipois, obs ", i)
    )
    expect_equal(
      log_lik_binegbin(i, mk(OLD, vint2 = V2)),
      log_lik_binegbin(i, mk(NEW, vint2 = V2)),
      label = paste0("binegbin_partialobs, obs ", i)
    )
  }
})

test_that("log_lik under the old names is non-degenerate", {
  # The values must be finite and must vary across draws and observations, so
  # that a fallback returning a constant would fail.
  ll <- vapply(seq_along(Y), function(i) log_lik_binegbin(i, mk(OLD)), numeric(3))
  expect_true(all(is.finite(ll)))
  expect_gt(length(unique(as.vector(ll))), 1)
})

test_that("posterior_epred is identical under old and new rate-dpar names", {
  expect_equal(
    posterior_epred_binegbin(mk(OLD)),
    posterior_epred_binegbin(mk(NEW))
  )
  expect_equal(
    posterior_epred_bipois(mk(OLD, five = FALSE)),
    posterior_epred_bipois(mk(NEW, five = FALSE))
  )
  # With vint2, posterior_epred_binegbin() resolves the rates through
  # .get_rate() and the exclusive dispersion through .SHAPEXTWO_NAMES, whose
  # last candidate is `shapex`. No bipois_partialobs() fit uses the old
  # spellings.
  expect_equal(
    posterior_epred_binegbin(mk(OLD, vint2 = V2)),
    posterior_epred_binegbin(mk(NEW, vint2 = V2))
  )
})

test_that("posterior_predict is identical under old and new rate-dpar names", {
  # The seed is reset before each call, so identical draws are expected.
  for (i in seq_along(Y)) {
    set.seed(42); a <- posterior_predict_binegbin(i, mk(OLD))
    set.seed(42); b <- posterior_predict_binegbin(i, mk(NEW))
    expect_identical(a, b, label = paste0("binegbin, obs ", i))

    set.seed(42); a <- posterior_predict_bipois(i, mk(OLD, five = FALSE))
    set.seed(42); b <- posterior_predict_bipois(i, mk(NEW, five = FALSE))
    expect_identical(a, b, label = paste0("bipois, obs ", i))

    set.seed(42); a <- posterior_predict_binegbin(i, mk(OLD, vint2 = V2))
    set.seed(42); b <- posterior_predict_binegbin(i, mk(NEW, vint2 = V2))
    expect_identical(a, b, label = paste0("binegbin_partialobs, obs ", i))
  }
})

test_that("the _y1_given_y2 and _y1_and_y2 workers are identical under old and new rate-dpar names", {
  # The workers in R/predict.R for the _y1_given_y2 functions and for
  # posterior_predict_y1_and_y2() read the rates through .get_rate().
  for (i in seq_along(Y)) {
    set.seed(42); a <- .impute_draw_binegbin(i, mk(OLD))
    set.seed(42); b <- .impute_draw_binegbin(i, mk(NEW))
    expect_identical(a, b, label = paste0("binegbin impute, obs ", i))

    set.seed(42); a <- .impute_draw_bipois(i, mk(OLD, five = FALSE))
    set.seed(42); b <- .impute_draw_bipois(i, mk(NEW, five = FALSE))
    expect_identical(a, b, label = paste0("bipois impute, obs ", i))

    set.seed(42); a <- .pair_draw_binegbin(i, mk(OLD))
    set.seed(42); b <- .pair_draw_binegbin(i, mk(NEW))
    expect_identical(a, b, label = paste0("binegbin pairs, obs ", i))

    set.seed(42); a <- .pair_draw_bipois(i, mk(OLD, five = FALSE))
    set.seed(42); b <- .pair_draw_bipois(i, mk(NEW, five = FALSE))
    expect_identical(a, b, label = paste0("bipois pairs, obs ", i))
  }
  expect_equal(.impute_mean_binegbin(mk(OLD)), .impute_mean_binegbin(mk(NEW)))
  expect_equal(.impute_mean_bipois(mk(OLD, five = FALSE)),
               .impute_mean_bipois(mk(NEW, five = FALSE)))
  expect_equal(.pair_mean(mk(OLD)), .pair_mean(mk(NEW)))

  # Under the old spelling and with the observation indicator, paired rows
  # return the recorded y1, with one column per observation.
  out <- .impute_y1_prep(mk(OLD, vint2 = V2), "binegbin",
                         type = "mean", keep_observed = TRUE)
  expect_equal(dim(out), c(3L, length(Y)))
  expect_equal(out[, V2 == 1L], matrix(Y[V2 == 1L], 3, sum(V2), byrow = TRUE))
})

test_that("a prep with neither rate-dpar spelling errors informatively", {
  bad <- make_synthetic_prep(
    dpars = list(mu = v$mu, lambdaX = v$one, lambdaY = v$two,
                 shapes = v$shapes, shapex = v$shapex),
    Y = Y, vint1 = V1
  )
  expect_error(log_lik_binegbin(1, bad), "lambdaone")
  expect_error(log_lik_binegbin(1, bad), "lambdaem")
  # The message lists the dpars that the prep declares.
  expect_error(log_lik_binegbin(1, bad), "lambdaX")
})

# ---------------------------------------------------------------------------
# Stored five-dpar binegbin fits
# ---------------------------------------------------------------------------
#
# Fits of binegbin made before 0.10.0 declare five dpars with a single
# `shapex`, and brms post-processes them with the current six-dpar methods
# because the family name is unchanged. Some stored fits declare
#
#   mu, lambdaem, lambdalb, shapes, shapex          (vint1 only, no vint2)
#
# which is the prep returned by mk(OLD). Three fallbacks must agree for such
# a fit: (1) .get_rate() resolves lambdaem and lambdalb; (2) .SHAPEXONE_NAMES
# and .SHAPEXTWO_NAMES both list `shapex` last; and (3) an absent vint2
# selects the paired branch. The tests above check (1); the tests below check
# (2), (3), and the combination of all three.

test_that("a stored five-dpar binegbin prep post-processes without error", {
  # The prep has no vint2 and declares `shapex` but not `shapexone`.
  stored <- mk(OLD)
  expect_null(stored$data$vint2)
  expect_true("shapex" %in% names(stored$dpars))
  expect_false("shapexone" %in% names(stored$dpars))

  ll <- vapply(seq_along(Y), function(i) log_lik_binegbin(i, stored), numeric(3))
  expect_true(all(is.finite(ll)))

  ep <- posterior_epred_binegbin(stored)
  expect_true(all(is.finite(ep)))
  expect_gt(length(unique(as.vector(ep))), 1)

  set.seed(7)
  pp <- vapply(seq_along(Y), function(i) posterior_predict_binegbin(i, stored), numeric(3))
  expect_true(all(pp >= 0))
})

test_that("the stored five-dpar answer equals the six-dpar answer with dispersions tied", {
  # A prep with a single `shapex` must give the same output as a prep with
  # shapexone == shapextwo == shapex, in every method.
  stored <- mk(OLD)
  tied <- make_synthetic_prep(
    dpars = list(mu = v$mu, lambdaone = v$one, lambdatwo = v$two,
                 shapes = v$shapes, shapexone = v$shapex, shapextwo = v$shapex),
    Y = Y, vint1 = V1
  )

  for (i in seq_along(Y)) {
    expect_equal(log_lik_binegbin(i, stored), log_lik_binegbin(i, tied),
                 label = paste0("log_lik, obs ", i))
    set.seed(42); a <- posterior_predict_binegbin(i, stored)
    set.seed(42); b <- posterior_predict_binegbin(i, tied)
    expect_identical(a, b, label = paste0("posterior_predict, obs ", i))
    # The marginal draw above reads shapexone only. The conditional draw and
    # the pair draw read both dispersions.
    set.seed(42); a <- .impute_draw_binegbin(i, stored)
    set.seed(42); b <- .impute_draw_binegbin(i, tied)
    expect_identical(a, b, label = paste0("impute draw, obs ", i))
    set.seed(42); a <- .pair_draw_binegbin(i, stored)
    set.seed(42); b <- .pair_draw_binegbin(i, tied)
    expect_identical(a, b, label = paste0("pair draw, obs ", i))
  }
  expect_equal(posterior_epred_binegbin(stored), posterior_epred_binegbin(tied))
  expect_equal(.impute_mean_binegbin(stored), .impute_mean_binegbin(tied))

  # Setting shapextwo to 0.05 changes log_lik and the conditional expectation,
  # so the checks above are not vacuous. At these rates, shapex = 5 is close to
  # the Poisson limit, and an 8-fold increase in shapextwo changes log_lik by
  # only about 4e-4.
  untied <- make_synthetic_prep(
    dpars = list(mu = v$mu, lambdaone = v$one, lambdatwo = v$two,
                 shapes = v$shapes, shapexone = v$shapex,
                 shapextwo = rep(0.05, length(v$shapex))),
    Y = Y, vint1 = V1
  )
  expect_gt(abs(log_lik_binegbin(1, untied)[1] - log_lik_binegbin(1, tied)[1]), 0.1)
  # The conditional expectation depends on shapextwo. The marginal mean
  # mu + lambdaone, returned here by posterior_epred_binegbin(), does not.
  expect_gt(max(abs(.impute_mean_binegbin(untied) - .impute_mean_binegbin(tied))), 0.1)
  expect_equal(posterior_epred_binegbin(untied), posterior_epred_binegbin(tied))
})
