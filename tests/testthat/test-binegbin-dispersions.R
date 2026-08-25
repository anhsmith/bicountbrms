# tests/testthat/test-binegbin-dispersions.R
#
# The two source-specific excess dispersions, `shapexone` and `shapextwo`.
#
# 0.8.0 split the single `shapex` into the pair, for the partially observed
# family only. 0.10.0 gives the pair to `binegbin()` as well -- which is the
# right way round, and the reason that release exists. Under partial
# observation `shapexone` is integrated out of the y2-only branch and so rests
# on the matched rows alone; under full pairing every row informs both, so that
# is the case in which they are easiest to identify. It was the case that
# lacked the capability.
#
# Recovery is therefore tested twice at the bottom of this file, once per
# constructor. Everything above it is shared: the likelihood is the same
# function either way, and the flag only chooses a branch.
#
# Four things need pinning, and they are separable:
#
#   1. The asymmetric likelihood is a proper pmf and keeps the marginal
#      identity -- i.e. freeing the dispersions did not break the
#      construction.
#   2. The STRUCTURE of the split: shapexone governs the y1 excess only, so it
#      must be absent from the y2-only branch, and swapping the two
#      dispersions must be equivalent to swapping the two rates.
#   3. Non-regression: with shapexone == shapextwo the six-dpar likelihood is
#      the pre-0.8.0 five-dpar one, and a five-dpar prep still resolves to both
#      excess dispersions without a shim.
#   4. Recovery of two genuinely different dispersions from simulated data.
#
# Everything here is self-contained: (1)-(3) need only R, (4) needs Stan.
# Agreement with the project-local family this release absorbs was a one-off
# migration check, not a test -- see the note above section 4.

# -----------------------------------------------------------------------
# 1. The asymmetric likelihood is a proper pmf
# -----------------------------------------------------------------------

# Summing a pmf over a TRUNCATED grid 0..K must land just below 1 and can
# never exceed it. Asserted one-sided for that reason: an upper bound of
# exactly 1 catches double-counting, which a symmetric tolerance would let
# through, while the lower bound absorbs the tail beyond K. The allowance has
# to be loose because a small dispersion is heavy-tailed -- at shapex = 0.4
# and lambda = 6 the mass past y = 150 is ~2e-5, and pushing K far enough to
# remove it costs K^2 evaluations of an inner sum.
expect_normalises <- function(mass, tol = 1e-4, label = NULL) {
  testthat::expect_lte(mass, 1 + 1e-9, label = label)
  testthat::expect_gt(mass, 1 - tol, label = label)
}

test_that("matched branch normalises to 1 with shapexone != shapextwo", {
  norm_check <- function(mu, lone, ltwo, ss, sx1, sx2, K = 150) {
    ys <- 0:K
    yg <- expand.grid(y1 = ys, y2 = ys)
    sum(exp(binegbin_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, ss, sx1, sx2)))
  }
  params <- list(
    c(mu = 5,  lone = 3, ltwo = 3, ss = 2,   sx1 = 0.5, sx2 = 8),
    c(mu = 10, lone = 2, ltwo = 6, ss = 1.5, sx1 = 20,  sx2 = 0.4),
    c(mu = 2,  lone = 4, ltwo = 1, ss = 5,   sx1 = 0.8, sx2 = 3)
  )
  for (p in params) {
    expect_normalises(
      norm_check(p[["mu"]], p[["lone"]], p[["ltwo"]], p[["ss"]], p[["sx1"]], p[["sx2"]]),
      label = paste(p, collapse = ",")
    )
  }
})

test_that("y2-only branch normalises to 1 over y2 with shapexone != shapextwo", {
  # One dimension only, so K can be pushed far enough for a tight bound.
  norm_check <- function(mu, ltwo, ss, sx1, sx2, K = 4000) {
    ys <- 0:K
    sum(exp(binegbin_lpmf_r(rep(0L, length(ys)), ys, y1_obs = 0L,
                                  mu, 3, ltwo, ss, sx1, sx2)))
  }
  expect_normalises(norm_check(5,  3, 2,   0.5, 8),   tol = 1e-8)
  expect_normalises(norm_check(10, 6, 1.5, 20,  0.4), tol = 1e-8)
})

test_that("moment identities hold with shapexone != shapextwo", {
  # The variance structure the README states for the asymmetric family. Each
  # margin picks up its OWN excess dispersion; the covariance picks up neither,
  # because only the shared component contributes to it. Computed by exact
  # summation over the joint pmf rather than by simulation, so this is an
  # identity check to machine precision, not a Monte Carlo one -- the pmf has
  # converged to 1 well inside the grid for these parameters (mass = 1 to
  # ~1e-16), which is what licenses the tight tolerance.
  #
  # A swapped shapexone/shapextwo would show up here as Vy1 and Vy2 trading
  # values, independently of the transposition check further down.
  nb_var <- function(m, phi) m + m^2 / phi

  check <- function(mu, l1, l2, ss, sx1, sx2, K = 150) {
    ys <- 0:K
    g  <- expand.grid(y1 = ys, y2 = ys)
    p  <- exp(binegbin_lpmf_r(g$y1, g$y2, mu, l1, l2, ss, sx1, sx2))
    Ey1 <- sum(p * g$y1); Ey2 <- sum(p * g$y2)
    lab <- paste(mu, l1, l2, ss, sx1, sx2, sep = ",")

    expect_equal(sum(p), 1, tolerance = 1e-10, label = paste("mass:", lab))
    expect_equal(Ey1, mu + l1, tolerance = 1e-8, label = paste("E[y1]:", lab))
    expect_equal(Ey2, mu + l2, tolerance = 1e-8, label = paste("E[y2]:", lab))
    expect_equal(sum(p * g$y1^2) - Ey1^2,
                 nb_var(mu, ss) + nb_var(l1, sx1),
                 tolerance = 1e-8, label = paste("Var(y1):", lab))
    expect_equal(sum(p * g$y2^2) - Ey2^2,
                 nb_var(mu, ss) + nb_var(l2, sx2),
                 tolerance = 1e-8, label = paste("Var(y2):", lab))
    # Cov depends on the shared component alone -- neither excess dispersion.
    expect_equal(sum(p * g$y1 * g$y2) - Ey1 * Ey2,
                 nb_var(mu, ss),
                 tolerance = 1e-8, label = paste("Cov:", lab))
    # Var(d) = the two excess variances summed; the shared count cancels.
    d <- g$y1 - g$y2
    expect_equal(sum(p * d^2) - sum(p * d)^2,
                 nb_var(l1, sx1) + nb_var(l2, sx2),
                 tolerance = 1e-8, label = paste("Var(d):", lab))
  }

  check(mu = 4, l1 = 3, l2 = 3, ss = 3, sx1 = 1.5, sx2 = 8)    # equal rates
  check(mu = 6, l1 = 2, l2 = 5, ss = 2, sx1 = 6,   sx2 = 1.2)  # asymmetry reversed
})

test_that("marginal identity holds under asymmetry", {
  # Integrating the matched branch over all y1 must reproduce the y2-only
  # branch. This is the identity that makes the unmatched rows coherent with
  # the matched ones, and it must survive the dispersions being freed.
  mu <- 6; lone <- 3; ltwo <- 4; ss <- 2; sx1 <- 0.6; sx2 <- 7
  K <- 400
  for (yl in c(0L, 1L, 3L, 7L, 15L)) {
    ys <- 0:K
    lp <- binegbin_lpmf_r(ys, rep(yl, length(ys)), mu, lone, ltwo, ss, sx1, sx2)
    marg <- { mx <- max(lp); mx + log(sum(exp(lp - mx))) }
    expect_equal(
      marg,
      binegbin_lpmf_r(0L, yl, y1_obs = 0L, mu, lone, ltwo, ss, sx1, sx2),
      tolerance = 1e-10, label = paste("y2 =", yl)
    )
  }
})

# -----------------------------------------------------------------------
# 2. Structure of the split
# -----------------------------------------------------------------------

test_that("shapexone does not enter the y2-only branch", {
  # The y1 margin is integrated out on y1_obs == 0 rows, taking its
  # dispersion with it. If shapexone leaked into that branch the two calls
  # below would differ -- and shapexone would be spuriously informed by the
  # unmatched rows, which is exactly the identifiability claim the docs make.
  base <- binegbin_lpmf_r(0L, 0:20, y1_obs = 0L, 6, 3, 4, 2, sx1 <- 0.5, 7)
  for (alt in c(0.01, 1, 100, 1e4)) {
    expect_identical(
      binegbin_lpmf_r(0L, 0:20, y1_obs = 0L, 6, 3, 4, 2, alt, 7),
      base, label = paste("shapexone =", alt)
    )
  }
  # ...and it DOES enter the matched branch, so the test above is not vacuous.
  m0 <- binegbin_lpmf_r(3L, 5L, 6, 3, 4, 2, 0.5, 7)
  m1 <- binegbin_lpmf_r(3L, 5L, 6, 3, 4, 2, 100, 7)
  expect_false(isTRUE(all.equal(m0, m1)))
})

test_that("swapping both rates and both dispersions transposes the joint", {
  # On the matched branch the two margins enter symmetrically, so
  # p(y1, y2 | lone, ltwo, sx1, sx2) == p(y2, y1 | ltwo, lone, sx2, sx1).
  # This pins that shapexone is bound to the y1 term and shapextwo to the y2
  # term, and not the reverse -- a swap the Stan signature could make
  # silently.
  ys <- 0:12
  yg <- expand.grid(y1 = ys, y2 = ys)
  a <- binegbin_lpmf_r(yg$y1, yg$y2, 6, 3, 4, 2, 0.6, 7)
  b <- binegbin_lpmf_r(yg$y2, yg$y1, 6, 4, 3, 2, 7, 0.6)
  expect_equal(a, b, tolerance = 1e-14)
})

# -----------------------------------------------------------------------
# 3. Non-regression on the symmetric special case
# -----------------------------------------------------------------------

test_that("shapexone == shapextwo reproduces the five-dpar likelihood exactly", {
  # The pre-0.8.0 family is this one under the constraint. Checked against
  # binegbin_lpmf_r (the single-dispersion reference, untouched by 0.8.0) on
  # the matched branch, and against the defaulted eight-argument call on
  # both.
  grid <- expand.grid(mu = c(0.5, 3, 12), lone = c(0.5, 2, 6),
                      ltwo = c(0.5, 2, 6), ss = c(0.8, 3), sx = c(0.8, 3))
  ys <- 0:15
  yg <- expand.grid(y1 = ys, y2 = ys)
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    expect_equal(
      binegbin_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, g$ss, g$sx, g$sx),
      binegbin_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, g$ss, g$sx),
      tolerance = 1e-14, label = paste(unlist(g), collapse = ",")
    )
  }
})

test_that("shapextwo defaults to shapexone (eight-argument calls are symmetric)", {
  ys <- 0:15
  yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
  expect_identical(
    binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, 6, 3, 4, 2, 1.5),
    binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, 6, 3, 4, 2, 1.5, 1.5)
  )
})

test_that("a five-dpar prep resolves `shapex` to BOTH excess dispersions", {
  # The stored-fit compatibility path, without needing a stored fit: a prep
  # carrying only the pre-0.8.0 `shapex` must behave exactly like one
  # carrying shapexone == shapextwo == that value. This is what lets
  # <= 0.7.0 fits post-process under 0.8.0 with no shim.
  Y <- c(5L, 9L, 0L, 14L); V1 <- c(4L, 7L, 2L, 11L); V2 <- c(1L, 0L, 1L, 0L)
  common <- list(mu = c(6, 6.5, 5.5), lambdaone = c(3, 3.2, 2.8),
                 lambdatwo = c(2, 2.1, 1.9), shapes = c(4, 4.5, 3.5))
  sx <- c(5, 5.5, 4.5)

  old <- make_synthetic_prep(c(common, list(shapex = sx)), Y, vint1 = V1, vint2 = V2)
  new <- make_synthetic_prep(c(common, list(shapexone = sx, shapextwo = sx)),
                             Y, vint1 = V1, vint2 = V2)
  for (i in seq_along(Y)) {
    expect_equal(log_lik_binegbin(i, old), log_lik_binegbin(i, new),
                 label = paste("obs", i))
    set.seed(11); a <- posterior_predict_binegbin(i, old)
    set.seed(11); b <- posterior_predict_binegbin(i, new)
    expect_identical(a, b, label = paste("posterior_predict, obs", i))
  }
})

test_that("a prep in the project-local ax spelling resolves to the package names", {
  # shapexem -> shapexone, shapexlb -> shapextwo. This is the mapping the
  # migration note documents, exercised in the direction that matters: an
  # already-stored asymmetric fit read by the packaged family.
  Y <- c(5L, 9L, 0L, 14L); V1 <- c(4L, 7L, 2L, 11L); V2 <- c(1L, 0L, 1L, 0L)
  common <- list(mu = c(6, 6.5, 5.5), shapes = c(4, 4.5, 3.5))
  one <- c(3, 3.2, 2.8); two <- c(2, 2.1, 1.9)
  sx1 <- c(0.5, 0.6, 0.55); sx2 <- c(7, 7.5, 6.5)

  ax  <- make_synthetic_prep(
    c(common, list(lambdaem = one, lambdalb = two, shapexem = sx1, shapexlb = sx2)),
    Y, vint1 = V1, vint2 = V2)
  pkg <- make_synthetic_prep(
    c(common, list(lambdaone = one, lambdatwo = two, shapexone = sx1, shapextwo = sx2)),
    Y, vint1 = V1, vint2 = V2)
  for (i in seq_along(Y)) {
    expect_equal(log_lik_binegbin(i, ax), log_lik_binegbin(i, pkg),
                 label = paste("obs", i))
  }
  # Not vacuous: the asymmetry is real in these numbers.
  flip <- make_synthetic_prep(
    c(common, list(lambdaone = one, lambdatwo = two, shapexone = sx2, shapextwo = sx1)),
    Y, vint1 = V1, vint2 = V2)
  expect_false(isTRUE(all.equal(log_lik_binegbin(1, pkg),
                                log_lik_binegbin(1, flip))))
})

test_that("a prep with no recognisable excess dispersion errors informatively", {
  bad <- make_synthetic_prep(
    list(mu = 6, lambdaone = 3, lambdatwo = 2, shapes = 4, shapeZ = 5),
    Y = 5L, vint1 = 4L, vint2 = 1L)
  expect_error(log_lik_binegbin(1, bad), "shapexone")
  expect_error(log_lik_binegbin(1, bad), "shapex")
  expect_error(log_lik_binegbin(1, bad), "shapeZ")
})

# -----------------------------------------------------------------------
# Stan: asymmetric grid against the R reference
# -----------------------------------------------------------------------

stan_code_asym <- paste0("functions {\n", binegbin_stan_funs, "}\nmodel {}\n")

stan_ready_asym <- FALSE
if (stan_tests_enabled() && requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = stan_code_asym)
      rstan::expose_stan_functions(sm)
    })
    stan_ready_asym <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan binegbin_lpmf matches the R reference with shapexone != shapextwo", {
  skip_if_not(stan_ready_asym, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(mu = c(0.2, 1, 5, 20), lone = c(0.2, 1, 5),
                      ltwo = c(0.2, 1, 5), ss = c(0.5, 2, 50),
                      sx1 = c(0.3, 2), sx2 = c(0.9, 30))
  check_one <- function(mu, lone, ltwo, ss, sx1, sx2) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
    stan_vals <- mapply(
      function(r, s, e) binegbin_lpmf(r, mu, lone, ltwo, ss, sx1, sx2, s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs,
                                    mu, lone, ltwo, ss, sx1, sx2)
    expect_true(all(is.finite(stan_vals)))
    max(abs(stan_vals - r_vals))
  }
  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, grid$ss,
                  grid$sx1, grid$sx2)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

# -----------------------------------------------------------------------
# 4. Simulation recovery with two genuinely different dispersions
# -----------------------------------------------------------------------

# NOT TESTED HERE: agreement with the project-local `binegbin_joint_ax()`
# family this release absorbs. That comparison was run once during the 0.8.0
# migration and is recorded in migration/family-unification.md; it is deliberately
# not a package test. It depends on fits stored outside this repository, so it
# could only ever skip on CI, on CRAN and on any other checkout -- a
# permanently-skipped test reads as coverage that is not there. Its subject is
# also frozen: those fits are fixed artefacts, so once the comparison passes it
# can only change if THIS package changes, and the package-side property is
# already pinned above without any external dependency (the Stan-vs-R
# asymmetric grid, the marginal identity, the transposition check, and
# shapexone's absence from the y2-only branch).

test_that("binegbin_partialobs recovers shapexone and shapextwo when they differ", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260803)

  n_vessel     <- 10L
  n_per_vessel <- 40L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int <- log(6)
  true_sd_vessel  <- 0.3
  true_lone       <- 6
  true_ltwo       <- 6
  true_shapes     <- 3
  true_shapexone  <- 0.7   # y1 excess: strongly overdispersed
  true_shapextwo  <- 6.0   # y2 excess: mildly overdispersed
  # An order of magnitude apart, so "both recovered" is a claim about the
  # split and not about a shared value being found twice.

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  z_mu   <- rnorm(n_vessel)
  mu_i   <- exp(true_log_mu_int + true_sd_vessel * z_mu[vessel])

  n_shared <- rnbinom(n, size = true_shapes,    mu = mu_i)
  n1       <- rnbinom(n, size = true_shapexone, mu = true_lone)
  n2       <- rnbinom(n, size = true_shapextwo, mu = true_ltwo)

  # shapexone is identified ONLY by the matched rows (it is integrated out of
  # the y2-only branch), so the design keeps a clear majority matched --
  # otherwise this test would measure the prior, not recovery.
  y1_obs <- rep(c(1L, 1L, 1L, 0L), length.out = n)

  dat <- data.frame(y1 = n_shared + n1, y2 = n_shared + n2,
                    y1_obs = y1_obs, vessel = factor(vessel))

  # brms leaves a custom family's non-mu dpars FLAT AND IMPROPER, and the
  # likelihood for a Negative-Binomial shape is flat towards large values (a
  # big shape is nearly Poisson, so the data cannot distinguish 30 from 300).
  # With one excess dispersion that is survivable; with two -- one of them,
  # shapexone, informed only by the matched rows -- it is not: run
  # unregularised this model produced 680 divergent transitions and a max Rhat
  # of 1.54, i.e. chains that never mixed. The README's "Set priors on the
  # dispersions" section says exactly this.
  #
  # normal(0, 1.5) on the log scale spans roughly [0.05, 20] at +/- 2 SD,
  # which contains all three true dispersions (0.7, 3, 6) without favouring
  # any of them -- it regularises the flat upper tail rather than supplying
  # the answer. Widening the SD would make it worse, not vaguer: these are log
  # links, so a wider normal is a heavier-tailed lognormal that puts MORE mass
  # in the flat region the sampler gets lost in.
  dispersion_priors <- c(
    brms::prior(normal(2, 1),   class = "Intercept"),
    brms::prior(normal(2, 1),   class = "b",         nlpar = "lamx"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapes"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapexone"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapextwo")
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2, y1_obs) ~ 1,
        mu ~ 1 + (1 | vessel),
        brms::nlf(lambdaone ~ lamx),
        brms::nlf(lambdatwo ~ lamx),
        lamx ~ 1, shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE
      ),
      family   = binegbin_partialobs(),
      stanvars = binegbin_partialobs_stanvars(),
      prior    = dispersion_priors,
      data     = dat,
      backend  = "rstan",
      chains   = 4, iter = 2000, warmup = 1000,
      seed     = 20260803, refresh = 0, init = 0.5,
      control  = list(adapt_delta = 0.99, max_treedepth = 12)
    )
  })

  draws <- as.data.frame(fit)

  # Smoke gates at a wide interval; see helper-coverage.R.
  expect_true(recovery_ok(draws, true_log_mu_int,      "b_Intercept"))
  expect_true(recovery_ok(draws, log(true_lone),       "b_lamx_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapes),     "b_shapes_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapexone),  "b_shapexone_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapextwo),  "b_shapextwo_Intercept"))

  # The substantive claim: the fit does not merely bracket each truth, it
  # resolves them as DIFFERENT. On the log scale the difference is
  # log(0.7) - log(6) = -2.15; a 95% interval excluding zero says the data
  # distinguished the two dispersions rather than the priors doing it.
  d  <- draws[["b_shapexone_Intercept"]] - draws[["b_shapextwo_Intercept"]]
  ci <- quantile(d, c(0.025, 0.975), names = FALSE)
  expect_lt(ci[[2]], 0)
  expect_true(
    log(true_shapexone) - log(true_shapextwo) >= ci[[1]] &&
      log(true_shapexone) - log(true_shapextwo) <= ci[[2]],
    label = paste0("true log-ratio ",
                   round(log(true_shapexone) - log(true_shapextwo), 3),
                   " outside [", round(ci[[1]], 3), ", ", round(ci[[2]], 3), "]")
  )

  # Dispatch must still resolve under six dpars.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))
  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))
  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# -----------------------------------------------------------------------
# 5. The same recovery under FULL pairing -- the capability 0.10.0 adds
# -----------------------------------------------------------------------
#
# The test above fits binegbin_partialobs() on a majority-matched design,
# because shapexone is identified only by the matched rows there. Under full
# pairing that constraint is gone: every row enters the matched branch and so
# informs both dispersions. This is the case the release exists to serve, and
# before 0.10.0 it could not be expressed at all -- binegbin() had one shapex,
# so there was nothing to recover separately.
#
# Same truths and the same priors as the partially observed fit, so the two are
# comparable. n is smaller (300 against 400) precisely because full pairing is
# the easier problem; if the design argument is right, the smaller fully paired
# design should still resolve the two.

test_that("binegbin recovers shapexone and shapextwo from a fully paired design", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260825)

  n_vessel     <- 10L
  n_per_vessel <- 30L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int <- log(6)
  true_sd_vessel  <- 0.3
  true_lone       <- 6
  true_ltwo       <- 6
  true_shapes     <- 3
  true_shapexone  <- 0.7   # y1 excess: strongly overdispersed
  true_shapextwo  <- 6.0   # y2 excess: mildly overdispersed

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  z_mu   <- rnorm(n_vessel)
  mu_i   <- exp(true_log_mu_int + true_sd_vessel * z_mu[vessel])

  n_shared <- rnbinom(n, size = true_shapes,    mu = mu_i)
  n1       <- rnbinom(n, size = true_shapexone, mu = true_lone)
  n2       <- rnbinom(n, size = true_shapextwo, mu = true_ltwo)

  # No y1_obs column anywhere. That is the point of the plain constructor: the
  # flag reaches Stan as a literal 1 and never enters the data.
  dat <- data.frame(y1 = n_shared + n1, y2 = n_shared + n2,
                    vessel = factor(vessel))

  dispersion_priors <- c(
    brms::prior(normal(2, 1),   class = "Intercept"),
    brms::prior(normal(2, 1),   class = "b",         nlpar = "lamx"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapes"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapexone"),
    brms::prior(normal(0, 1.5), class = "Intercept", dpar  = "shapextwo")
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2) ~ 1,
        mu ~ 1 + (1 | vessel),
        brms::nlf(lambdaone ~ lamx),
        brms::nlf(lambdatwo ~ lamx),
        lamx ~ 1, shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE
      ),
      family   = binegbin(),
      stanvars = binegbin_stanvars(),
      prior    = dispersion_priors,
      data     = dat,
      backend  = "rstan",
      chains   = 4, iter = 2000, warmup = 1000,
      seed     = 20260825, refresh = 0, init = 0.5,
      control  = list(adapt_delta = 0.99, max_treedepth = 12)
    )
  })

  # The generated model must carry six dpars and no vint2.
  sd_fit <- brms::standata(fit)
  expect_true("vint1" %in% names(sd_fit))
  expect_false("vint2" %in% names(sd_fit))

  draws <- as.data.frame(fit)

  expect_true(recovery_ok(draws, true_log_mu_int,     "b_Intercept"))
  expect_true(recovery_ok(draws, log(true_lone),      "b_lamx_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapes),    "b_shapes_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapexone), "b_shapexone_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapextwo), "b_shapextwo_Intercept"))

  # The substantive claim, as in the partially observed fit: the two are
  # resolved as DIFFERENT, not bracketed twice around a shared value.
  d  <- draws[["b_shapexone_Intercept"]] - draws[["b_shapextwo_Intercept"]]
  ci <- quantile(d, c(0.025, 0.975), names = FALSE)
  expect_lt(ci[[2]], 0)
  expect_true(
    log(true_shapexone) - log(true_shapextwo) >= ci[[1]] &&
      log(true_shapexone) - log(true_shapextwo) <= ci[[2]],
    label = paste0("true log-ratio ",
                   round(log(true_shapexone) - log(true_shapextwo), 3),
                   " outside [", round(ci[[1]], 3), ", ", round(ci[[2]], 3), "]")
  )

  # Dispatch resolves to the unified methods from a fit with no vint2.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))
  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))
  ep <- brms::posterior_epred(fit)
  expect_equal(dim(ep)[2], n)
  expect_true(all(is.finite(ep)))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))
  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# NOT COVERED, deliberately. The design argument above says a fully paired
# design identifies shapexone MORE sharply than a partially paired one of the
# same size. That is a directional claim about two posterior SDs, and testing
# it would mean fitting both at matched n and comparing -- two Stan fits to
# assert an inequality whose margin is not known in advance. The claim is
# documented in ?binegbin_partialobs and in the README; what is pinned here is
# the weaker and more useful fact that both designs recover both dispersions.
