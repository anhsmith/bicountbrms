# tests/testthat/test-binegbin-dispersions.R
#
# The dispersions `shapexone` and `shapextwo` of the two source-exclusive
# components. Sections 1-3 check, in R alone, that the likelihood with unequal
# dispersions is a proper pmf, that `shapexone` is absent from the unpaired
# branch, and that equal dispersions reproduce the single-dispersion
# likelihood. Sections 4 and 5 fit models in Stan and check recovery under
# partial and full pairing.

# -----------------------------------------------------------------------
# 1. The asymmetric likelihood is a proper pmf
# -----------------------------------------------------------------------

# A pmf summed over the truncated grid 0..K cannot exceed 1, so the upper
# bound is 1 (which detects double-counting) and the lower bound allows for
# the mass beyond K. At a dispersion of 0.4 and a rate of 6, the mass beyond
# y = 150 is about 2e-5, so the default lower tolerance is 1e-4.
expect_normalises <- function(mass, tol = 1e-4, label = NULL) {
  testthat::expect_lte(mass, 1 + 1e-9, label = label)
  testthat::expect_gt(mass, 1 - tol, label = label)
}

test_that("paired branch normalises to 1 with shapexone != shapextwo", {
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

test_that("unpaired branch normalises to 1 over y2 with shapexone != shapextwo", {
  # The sum runs over y2 alone, so K = 4000 is affordable and the bound is
  # tight.
  norm_check <- function(mu, ltwo, ss, sx1, sx2, K = 4000) {
    ys <- 0:K
    sum(exp(binegbin_lpmf_r(rep(0L, length(ys)), ys, y1_obs = 0L,
                                  mu, 3, ltwo, ss, sx1, sx2)))
  }
  expect_normalises(norm_check(5,  3, 2,   0.5, 8),   tol = 1e-8)
  expect_normalises(norm_check(10, 6, 1.5, 20,  0.4), tol = 1e-8)
})

test_that("moment identities hold with shapexone != shapextwo", {
  # Var(y1) depends on shapexone, Var(y2) on shapextwo, and Cov(y1, y2) on
  # neither exclusive dispersion. The moments are computed by exact summation
  # over the joint pmf, whose mass on this grid equals 1 to within 1e-16, so
  # the tolerances are close to machine precision.
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
    # The covariance depends on the shared component alone.
    expect_equal(sum(p * g$y1 * g$y2) - Ey1 * Ey2,
                 nb_var(mu, ss),
                 tolerance = 1e-8, label = paste("Cov:", lab))
    # Var(d) is the sum of the variances of the two exclusive components,
    # because the shared count cancels from d = y1 - y2.
    d <- g$y1 - g$y2
    expect_equal(sum(p * d^2) - sum(p * d)^2,
                 nb_var(l1, sx1) + nb_var(l2, sx2),
                 tolerance = 1e-8, label = paste("Var(d):", lab))
  }

  check(mu = 4, l1 = 3, l2 = 3, ss = 3, sx1 = 1.5, sx2 = 8)    # equal rates
  check(mu = 6, l1 = 2, l2 = 5, ss = 2, sx1 = 6,   sx2 = 1.2)  # asymmetry reversed
})

test_that("marginal identity holds under asymmetry", {
  # The paired branch summed over y1 must equal the unpaired branch, with
  # unequal dispersions.
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

test_that("shapexone does not enter the unpaired branch", {
  # On rows with y1_obs == 0, the y1 margin is summed out, so the unpaired
  # branch must not depend on shapexone.
  base <- binegbin_lpmf_r(0L, 0:20, y1_obs = 0L, 6, 3, 4, 2, sx1 <- 0.5, 7)
  for (alt in c(0.01, 1, 100, 1e4)) {
    expect_identical(
      binegbin_lpmf_r(0L, 0:20, y1_obs = 0L, 6, 3, 4, 2, alt, 7),
      base, label = paste("shapexone =", alt)
    )
  }
  # The paired branch depends on shapexone, so the check above is not vacuous.
  m0 <- binegbin_lpmf_r(3L, 5L, 6, 3, 4, 2, 0.5, 7)
  m1 <- binegbin_lpmf_r(3L, 5L, 6, 3, 4, 2, 100, 7)
  expect_false(isTRUE(all.equal(m0, m1)))
})

test_that("swapping both rates and both dispersions transposes the joint", {
  # On the paired branch,
  # p(y1, y2 | lone, ltwo, sx1, sx2) == p(y2, y1 | ltwo, lone, sx2, sx1).
  # The identity fails if shapexone and shapextwo are bound to the wrong
  # exclusive components.
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
  # On the paired branch, a call with shapexone == shapextwo must equal a call
  # that supplies shapexone alone.
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

test_that("a five-dpar prep resolves `shapex` to BOTH exclusive-component dispersions", {
  # Stored fits made before 0.8.0 declare a single `shapex`. A prep with
  # `shapex` must give the same log_lik and posterior_predict as a prep with
  # shapexone == shapextwo == shapex.
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

test_that("a prep with the dpar names lambdaem, lambdalb, shapexem, and shapexlb resolves to the current names", {
  # Some stored fits declare the dpar names lambdaem, lambdalb, shapexem, and
  # shapexlb. These must resolve to lambdaone, lambdatwo, shapexone, and
  # shapextwo, in that order.
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
  # Swapping the two dispersions changes log_lik, so the check above is not
  # vacuous.
  flip <- make_synthetic_prep(
    c(common, list(lambdaone = one, lambdatwo = two, shapexone = sx2, shapextwo = sx1)),
    Y, vint1 = V1, vint2 = V2)
  expect_false(isTRUE(all.equal(log_lik_binegbin(1, pkg),
                                log_lik_binegbin(1, flip))))
})

test_that("a prep with no recognisable exclusive-component dispersion errors informatively", {
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
# 4. Recovery of two different dispersions under partial pairing
# -----------------------------------------------------------------------

test_that("binegbin_partialobs recovers shapexone and shapextwo when they differ", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260803)

  n_group     <- 10L
  n_per_group <- 40L
  n            <- n_group * n_per_group

  true_log_mu_int <- log(6)
  true_sd_group  <- 0.3
  true_lone       <- 6
  true_ltwo       <- 6
  true_shapes     <- 3
  true_shapexone  <- 0.7   # y1 exclusive component: strongly overdispersed
  true_shapextwo  <- 6.0   # y2 exclusive component: mildly overdispersed
  # The two dispersions differ by an order of magnitude, so recovery of both
  # cannot come from one shared value.

  group <- rep(seq_len(n_group), each = n_per_group)
  z_mu   <- rnorm(n_group)
  mu_i   <- exp(true_log_mu_int + true_sd_group * z_mu[group])

  n_shared <- rnbinom(n, size = true_shapes,    mu = mu_i)
  n1       <- rnbinom(n, size = true_shapexone, mu = true_lone)
  n2       <- rnbinom(n, size = true_shapextwo, mu = true_ltwo)

  # Only the paired rows inform shapexone, so three rows in four are paired.
  y1_obs <- rep(c(1L, 1L, 1L, 0L), length.out = n)

  dat <- data.frame(y1 = n_shared + n1, y2 = n_shared + n2,
                    y1_obs = y1_obs, group = factor(group))

  # brms gives the non-mu dpars of a custom family flat improper priors, and
  # the likelihood of a negative binomial dispersion is flat towards large
  # values. Without priors, this model gave 680 divergent transitions and a
  # maximum Rhat of 1.54. A normal(0, 1.5) prior on the log scale spans about
  # 0.05 to 20 at +/- 2 SD, which includes the three true dispersions (0.7, 3,
  # and 6). A wider SD on the log scale would place more mass at large
  # dispersions, where the likelihood is flat.
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
        mu ~ 1 + (1 | group),
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

  # The 95% interval for log(shapexone) - log(shapextwo) must exclude zero and
  # contain the true value, log(0.7) - log(6) = -2.15. The priors on the two
  # dispersions are identical, so the separation comes from the data.
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

  # Post-processing methods dispatch on the fitted model with six dpars.
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
# 5. Recovery of two different dispersions under full pairing
# -----------------------------------------------------------------------
#
# Under full pairing, every row informs both shapexone and shapextwo. The true
# values and priors are those of section 4, with 300 rows rather than 400.

test_that("binegbin recovers shapexone and shapextwo from a fully paired design", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260825)

  n_group     <- 10L
  n_per_group <- 30L
  n            <- n_group * n_per_group

  true_log_mu_int <- log(6)
  true_sd_group  <- 0.3
  true_lone       <- 6
  true_ltwo       <- 6
  true_shapes     <- 3
  true_shapexone  <- 0.7   # y1 exclusive component: strongly overdispersed
  true_shapextwo  <- 6.0   # y2 exclusive component: mildly overdispersed

  group <- rep(seq_len(n_group), each = n_per_group)
  z_mu   <- rnorm(n_group)
  mu_i   <- exp(true_log_mu_int + true_sd_group * z_mu[group])

  n_shared <- rnbinom(n, size = true_shapes,    mu = mu_i)
  n1       <- rnbinom(n, size = true_shapexone, mu = true_lone)
  n2       <- rnbinom(n, size = true_shapextwo, mu = true_ltwo)

  # binegbin() passes the observation indicator to Stan as the literal 1, so
  # the data contain no y1_obs column.
  dat <- data.frame(y1 = n_shared + n1, y2 = n_shared + n2,
                    group = factor(group))

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
        mu ~ 1 + (1 | group),
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

  # The Stan data of a binegbin() fit contain vint1 and no vint2.
  sd_fit <- brms::standata(fit)
  expect_true("vint1" %in% names(sd_fit))
  expect_false("vint2" %in% names(sd_fit))

  draws <- as.data.frame(fit)

  expect_true(recovery_ok(draws, true_log_mu_int,     "b_Intercept"))
  expect_true(recovery_ok(draws, log(true_lone),      "b_lamx_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapes),    "b_shapes_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapexone), "b_shapexone_Intercept"))
  expect_true(recovery_ok(draws, log(true_shapextwo), "b_shapextwo_Intercept"))

  # As in section 4, the 95% interval for the log ratio of the dispersions
  # must exclude zero and contain the true value.
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

  # Post-processing methods dispatch on a fit with no vint2.
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
