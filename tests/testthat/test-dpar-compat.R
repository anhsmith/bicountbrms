# Regression tests for the pre-0.7.0 dpar-name fallback (.get_rate).
#
# 0.7.0 renamed the joint families' excess-rate dpars lambdaem/lambdalb to
# lambdaone/lambdatwo. A brmsfit stores its OWN family object, so
# prepare_predictions() on a fit made before the rename hands post-processing a
# prep whose dpars carry the old names. Every rate read in bipois/bipois_partialobs/
# binegbin/binegbin_partialobs therefore goes through .get_rate(), which resolves the name
# against the fit rather than assuming the current spelling.
#
# What these tests pin: identical inputs under either spelling produce
# identical output (so old fits are not merely non-erroring but numerically
# unchanged), and a prep carrying neither spelling fails with a message that
# names the problem instead of a bare `!is.null(x) is not TRUE` from
# brms::get_dpar().
#
# These use make_synthetic_prep() (helper-prep.R) rather than a real fit, so
# they run without Stan. The real-fit path -- where brms also resolves each
# dpar's LINK by name via prep$family -- is exercised separately against an
# actual pre-rename fit.

# Length-3 vectors => 3 posterior draws; 4 observations.
v <- list(mu = c(6, 6.5, 5.5), one = c(3, 3.2, 2.8), two = c(2, 2.1, 1.9),
          shapes = c(4, 4.5, 3.5), shapex = c(5, 5.5, 4.5))
Y  <- c(5L, 9L, 0L, 14L)
V1 <- c(4L, 7L, 2L, 11L)
V2 <- c(1L, 0L, 1L, 0L)

# Same numbers, two spellings. `nm` picks which pair of rate names to use.
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
  # Guards against the fallback silently returning a constant (e.g. if the rate
  # resolved to something wrong but finite): the values must be finite and must
  # vary across draws and observations.
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
  # binegbin_partialobs gained a posterior_epred at 0.9.0, so a pre-existing fit is
  # post-processed by code it never ran under. It resolves both rates through
  # .get_rate() and its excess dispersion through .SHAPEXTWO_NAMES, whose last
  # candidate is the five-dpar `shapex` supplied by mk() -- the same fallback
  # log_lik and posterior_predict already rely on. bipois_partialobs is not tested
  # here: it is new at 0.9.0, so no fit can carry the old spellings.
  expect_equal(
    posterior_epred_binegbin(mk(OLD, vint2 = V2)),
    posterior_epred_binegbin(mk(NEW, vint2 = V2))
  )
})

test_that("posterior_predict is identical under old and new rate-dpar names", {
  # Stochastic, so the seed is reset before each call rather than compared
  # distributionally.
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

test_that("a prep with neither rate-dpar spelling errors informatively", {
  bad <- make_synthetic_prep(
    dpars = list(mu = v$mu, lambdaX = v$one, lambdaY = v$two,
                 shapes = v$shapes, shapex = v$shapex),
    Y = Y, vint1 = V1
  )
  expect_error(log_lik_binegbin(1, bad), "lambdaone")
  expect_error(log_lik_binegbin(1, bad), "lambdaem")
  # The message should list what the fit actually has, so the user can see why.
  expect_error(log_lik_binegbin(1, bad), "lambdaX")
})

# ---------------------------------------------------------------------------
# The shape a stored plain-binegbin fit actually has
# ---------------------------------------------------------------------------
#
# 0.10.0 keeps the name `binegbin` and changes what it means: five dpars with a
# single `shapex` become six with shapexone/shapextwo. A fit stored under the
# old meaning therefore lands on the NEW post-processing methods no matter what
# anyone pins, because brms resolves log_lik_<name> off the attached search
# path using the name the fit stored.
#
# Such fits exist. Scanning tnc001-belize-em/fits at 0.9.1 found 103 with
# family$name == "binegbin", and every one of them carries
#
#   mu, lambdaem, lambdalb, shapes, shapex          (vint1 only, no vint2)
#
# -- pre-0.7.0 RATE names as well as the single dispersion. Three independent
# fallbacks have to line up for those to keep working:
#
#   1. .get_rate() resolves lambdaem/lambdalb to the one/two positions;
#   2. .SHAPEXONE_NAMES and .SHAPEXTWO_NAMES both list `shapex` last, so both
#      per-margin dispersions resolve to the one the fit has;
#   3. an absent vint2 selects the matched branch.
#
# mk(OLD) is that prep exactly. The blocks above already pin (1); these pin
# (2) and (3), and pin the composition -- which is the part that holds only
# because three separate mechanisms happen to agree.

test_that("a stored five-dpar binegbin prep post-processes without error", {
  # Pins the shape itself, independently of what it is compared against.
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
  # The composition test. A single `shapex` must mean exactly what supplying
  # shapexone == shapextwo == shapex means, in all three methods. If either
  # .SHAPEX*_NAMES vector stopped falling through to `shapex`, or fell through
  # to a different dpar, this is what would catch it.
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
  }
  expect_equal(posterior_epred_binegbin(stored), posterior_epred_binegbin(tied))

  # Not vacuous: untying the two dispersions changes both. shapextwo is taken
  # two orders below shapex rather than a small multiple above it -- at these
  # rates the likelihood is nearly flat in the direction of LESS
  # overdispersion (shapex = 5 is already close to the Poisson limit for a
  # mean of 2), so an 8x increase moves log_lik by only ~4e-4 and would make
  # this assertion a coin toss.
  untied <- make_synthetic_prep(
    dpars = list(mu = v$mu, lambdaone = v$one, lambdatwo = v$two,
                 shapes = v$shapes, shapexone = v$shapex,
                 shapextwo = rep(0.05, length(v$shapex))),
    Y = Y, vint1 = V1
  )
  expect_gt(abs(log_lik_binegbin(1, untied)[1] - log_lik_binegbin(1, tied)[1]), 0.1)
  expect_gt(max(abs(posterior_epred_binegbin(untied) - posterior_epred_binegbin(tied))), 0.1)
})
