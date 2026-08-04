# Regression tests for the pre-0.7.0 dpar-name fallback (.get_rate).
#
# 0.7.0 renamed the joint families' excess-rate dpars lambdaem/lambdalb to
# lambdaone/lambdatwo. A brmsfit stores its OWN family object, so
# prepare_predictions() on a fit made before the rename hands post-processing a
# prep whose dpars carry the old names. Every rate read in bipois/bipois_cens/
# binegbin/binegbin_cens therefore goes through .get_rate(), which resolves the name
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
      log_lik_binegbin_cens(i, mk(OLD, vint2 = V2)),
      log_lik_binegbin_cens(i, mk(NEW, vint2 = V2)),
      label = paste0("binegbin_cens, obs ", i)
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
  # binegbin_cens gained a posterior_epred at 0.9.0, so a pre-existing fit is
  # post-processed by code it never ran under. It resolves both rates through
  # .get_rate() and its excess dispersion through .SHAPEXTWO_NAMES, whose last
  # candidate is the five-dpar `shapex` supplied by mk() -- the same fallback
  # log_lik and posterior_predict already rely on. bipois_cens is not tested
  # here: it is new at 0.9.0, so no fit can carry the old spellings.
  expect_equal(
    posterior_epred_binegbin_cens(mk(OLD, vint2 = V2)),
    posterior_epred_binegbin_cens(mk(NEW, vint2 = V2))
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

    set.seed(42); a <- posterior_predict_binegbin_cens(i, mk(OLD, vint2 = V2))
    set.seed(42); b <- posterior_predict_binegbin_cens(i, mk(NEW, vint2 = V2))
    expect_identical(a, b, label = paste0("binegbin_cens, obs ", i))
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
