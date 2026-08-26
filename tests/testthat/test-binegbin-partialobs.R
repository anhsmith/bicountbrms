# tests/testthat/test-binegbin-partialobs.R

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("matched (y1_obs==1) joint PMF normalises to 1 across parameter sets", {
  norm_check <- function(mu, lone, ltwo, ss, sx, K = 120) {
    ys <- 0:K
    yg <- expand.grid(y1 = ys, y2 = ys)
    lp <- binegbin_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, ss, sx)
    sum(exp(lp))
  }
  params <- list(
    c(mu = 5,  lone = 3, ltwo = 3, ss = 2,   sx = 1),
    c(mu = 10, lone = 2, ltwo = 2, ss = 1.5, sx = 2),
    c(mu = 2,  lone = 4, ltwo = 4, ss = 5,   sx = 0.8)
  )
  for (p in params) {
    expect_equal(norm_check(p[["mu"]], p[["lone"]], p[["ltwo"]], p[["ss"]], p[["sx"]]),
                 1, tolerance = 1e-6, label = paste(p, collapse = ","))
  }
})

test_that("LB-only (y1_obs==0) branch normalises to 1 over y2", {
  # The y1_obs==0 branch is a 1-D pmf in y2 (the y1 margin is integrated
  # out), so it must sum to 1 over y2 alone.
  norm_check <- function(mu, lone, ltwo, ss, sx, K = 200) {
    ys <- 0:K
    lp <- binegbin_lpmf_r(rep(0L, length(ys)), ys, y1_obs = 0L, mu, lone, ltwo, ss, sx)
    sum(exp(lp))
  }
  params <- list(
    c(mu = 5,  lone = 3, ltwo = 3, ss = 2,   sx = 1),
    c(mu = 10, lone = 2, ltwo = 6, ss = 1.5, sx = 2),
    c(mu = 2,  lone = 4, ltwo = 4, ss = 5,   sx = 0.8)
  )
  for (p in params) {
    expect_equal(norm_check(p[["mu"]], p[["lone"]], p[["ltwo"]], p[["ss"]], p[["sx"]]),
                 1, tolerance = 1e-6, label = paste(p, collapse = ","))
  }
})

test_that("marginal identity: sum over y1 of the matched branch == LB-only branch", {
  # Integrating the matched (y1_obs==1) joint over all y1 must reproduce the
  # LB-only (y1_obs==0) branch value for that y2, exactly (both are the
  # y1-integrated marginal of the same bivariate model).
  mu <- 6; lone <- 3; ltwo <- 4; ss <- 2; sx <- 1.5
  K  <- 300
  y2_vals <- c(0L, 1L, 3L, 7L, 15L)
  for (yl in y2_vals) {
    ys <- 0:K
    lp_joint <- binegbin_lpmf_r(ys, rep(yl, length(ys)), mu, lone, ltwo, ss, sx)
    marg_from_joint <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
    lb_branch <- binegbin_lpmf_r(0L, yl, y1_obs = 0L, mu, lone, ltwo, ss, sx)
    expect_equal(marg_from_joint, lb_branch, tolerance = 1e-10,
                 label = paste("y2 =", yl))
  }
})

# Up to 0.9.1 a block here checked that binegbin_partialobs's matched branch was
# byte-for-byte the binegbin lpmf, because they were two separate
# implementations that could drift. Since 0.10.0 there is one lpmf and the
# identity is structural, so that check would compare a function to itself.
# What is still worth pinning is which branch an UNFLAGGED call selects --
# that is a default, not a structural fact, and getting it wrong would silently
# score every fully paired row against the marginal.

test_that("an unflagged call is the matched branch, not the marginal", {
  grid <- expand.grid(
    mu   = c(0.5, 3, 12),
    lone = c(0.5, 2, 6),
    ltwo = c(0.5, 2, 6),
    ss   = c(0.8, 3),
    sx   = c(0.8, 3)
  )
  ys <- 0:20
  yg <- expand.grid(y1 = ys, y2 = ys)
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    default  <- binegbin_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, g$ss, g$sx)
    explicit <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = 1L,
                                g$mu, g$lone, g$ltwo, g$ss, g$sx)
    expect_equal(default, explicit, tolerance = 1e-14,
                 label = paste(unlist(g), collapse = ","))
  }
  # Not vacuous: the other branch is a different number on almost every row.
  g <- grid[1, ]
  marg <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = 0L,
                          g$mu, g$lone, g$ltwo, g$ss, g$sx)
  base <- binegbin_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, g$ss, g$sx)
  expect_gt(mean(abs(marg - base) > 1e-8), 0.9)
})

test_that("posterior_predict draws reproduce the joint/marginal conditional y1 | y2", {
  # posterior_predict_binegbin samples N_shared | y2 then adds a fresh
  # N10. Its distribution must equal P(y1 | y2) = joint(y1, y2) /
  # marginal(y2) -- the exact discrete conditional. Checked by Monte Carlo.
  set.seed(20260706)
  mu <- 5; lone <- 3; ltwo <- 4; ss <- 2; sx <- 1.5
  y2 <- 6L
  ndraws <- 2e5

  prep <- make_synthetic_prep(
    dpars = list(
      mu       = rep(mu,  ndraws),
      lambdaone = rep(lone, ndraws),
      lambdatwo = rep(ltwo, ndraws),
      shapes   = rep(ss,  ndraws),
      shapex   = rep(sx,  ndraws)
    ),
    Y     = 0L,          # response value is unused by posterior_predict
    vint1 = y2,
    vint2 = 1L
  )
  draws <- posterior_predict_binegbin(1, prep)
  expect_length(draws, ndraws)

  # Analytic conditional P(y1 = x | y2) = joint(x, y2) / marginal(y2).
  K <- 80
  xs <- 0:K
  lp_joint <- binegbin_lpmf_r(xs, rep(y2, length(xs)), mu, lone, ltwo, ss, sx)
  lp_marg  <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
  p_cond   <- exp(lp_joint - lp_marg)

  emp <- tabulate(draws + 1L, nbins = K + 1L) / ndraws
  # Compare on the well-populated support; MC error ~ 1/sqrt(ndraws).
  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

# One lpmf, so one function to compile. Up to 0.9.1 this compiled binegbin and
# binegbin_partialobs together, because the matched branch of the second was a second
# copy of the first and the two could drift; the Stan-level equivalence check
# existed to stop that. Since 0.10.0 both constructors call this one function,
# so the equivalence holds by construction and there is nothing to compare.
# What remains is the grid check of Stan against the independent R reference,
# on both branches.
stan_code <- paste0("functions {\n", binegbin_stan_funs, "}\nmodel {}\n")

stan_ready <- FALSE
if (stan_tests_enabled() && requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = stan_code)
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan binegbin_lpmf matches R brute-force reference (both branches)", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu  = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5),
    ss  = c(0.5, 2, 50),
    sx  = c(0.5, 2, 50)
  )

  check_one <- function(mu, lone, ltwo, ss, sx) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
    # Symmetric case: the single `sx` is passed for BOTH shapexone and
    # shapextwo, which is the constraint the pre-0.8.0 five-dpar family
    # imposed. The asymmetric grid is in test-binegbin-dispersions.R.
    stan_vals <- mapply(
      function(r, s, e) binegbin_lpmf(r, mu, lone, ltwo, ss, sx, sx, s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, mu, lone, ltwo, ss, sx)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, grid$ss, grid$sx)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("the Stan matched branch routes the two dispersions the way R does", {
  # The Stan-vs-Stan equivalence check this replaces compared binegbin_lpmf
  # with the 0.9.1 binegbin_partialobs_lpmf; there is now one function, so it would compare a
  # function to itself. The hazard it guarded against -- the two
  # implementations disagreeing -- is gone, but a sharper one remains: the Stan
  # signature and the R reference's argument order can be reordered
  # independently, silently swapping which dispersion governs which component.
  # Checked here with shapexone far from shapextwo, which the grid above does
  # not do.
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5),
    ss   = c(0.5, 2, 50)
  )
  SX1 <- 0.6; SX2 <- 7      # an order of magnitude apart, deliberately

  check_one <- function(mu, lone, ltwo, ss) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(
      function(r, s) binegbin_lpmf(r, mu, lone, ltwo, ss, SX1, SX2, s, 1L),
      yg$y1, yg$y2)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, ss, SX1, SX2)
    max(abs(stan_vals - r_vals))
  }
  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, grid$ss)
  expect_true(max(diffs) < 1e-12, label = paste("max diff =", max(diffs)))

  # Not vacuous: handing Stan the two dispersions the other way round is a
  # materially different number, so the agreement above is a real constraint
  # on the argument order rather than an artefact of the two being close.
  swapped <- mapply(function(r, s) binegbin_lpmf(r, 5, 2, 3, 2, SX2, SX1, s, 1L),
                    c(9L, 4L, 12L), c(6L, 7L, 3L))
  straight <- binegbin_lpmf_r(c(9L, 4L, 12L), c(6L, 7L, 3L), 5, 2, 3, 2, SX1, SX2)
  expect_gt(max(abs(swapped - straight)), 0.05)
})

test_that("Stan binegbin_lpmf is numerically stable at extreme rates and shapes", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  edge_cases <- list(
    c(mu = 1e-6, lone = 1e-6, ltwo = 1e-6, ss = 0.2, sx = 0.2),
    c(mu = 1e-6, lone = 50,   ltwo = 50,   ss = 100, sx = 100),
    c(mu = 200,  lone = 1e-6, ltwo = 1e-6, ss = 0.2, sx = 0.2),
    c(mu = 200,  lone = 200,  ltwo = 200,  ss = 50,  sx = 50)
  )
  for (ec in edge_cases) {
    ys <- 0:5
    yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
    stan_vals <- mapply(
      function(r, s, e) binegbin_lpmf(r, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], ec[["ss"]], ec[["sx"]], ec[["sx"]], s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs,
                                    ec[["mu"]], ec[["lone"]], ec[["ltwo"]], ec[["ss"]], ec[["sx"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# brms end-to-end: dispatch (loo + posterior_predict) and recovery
# -----------------------------------------------------------------------

test_that("binegbin_partialobs fits the SYMMETRIC model via nlf, dispatches, and recovers params", {
  # Since 0.8.0 the two excess dispersions are separate dpars, so the
  # symmetric model -- one dispersion for both excess components, which is
  # what this simulation generates -- is expressed as a formula constraint:
  # both dpars routed through a single non-linear parameter `shapexx`. This
  # is the migration idiom documented in migration/family-unification.md, so the
  # test doubles as a check that it fits and recovers.
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260706)

  n_vessel     <- 8L
  n_per_vessel <- 30L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int <- log(8)
  true_sd_vessel  <- 0.3
  true_lone        <- 3           # shared excess rate (lambdaone = lambdatwo)
  true_shapes     <- 2
  true_shapex     <- 1.5

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  z_mu   <- rnorm(n_vessel)
  mu_i   <- exp(true_log_mu_int + true_sd_vessel * z_mu[vessel])

  n_shared <- rnbinom(n, size = true_shapes, mu = mu_i)
  n1      <- rnbinom(n, size = true_shapex, mu = true_lone)
  n2      <- rnbinom(n, size = true_shapex, mu = true_lone)

  # Half the rows are LB-only (y1 unobserved) -- the partial observation the family
  # exists to handle.
  y1_obs <- rep(c(1L, 0L), length.out = n)

  dat <- data.frame(
    y1   = n_shared + n1,
    y2   = n_shared + n2,
    y1_obs = y1_obs,
    vessel = factor(vessel)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2, y1_obs) ~ 1,
        mu ~ 1 + (1 | vessel),
        brms::nlf(lambdaone ~ lamx),
        brms::nlf(lambdatwo ~ lamx),
        brms::nlf(shapexone ~ shapexx),
        brms::nlf(shapextwo ~ shapexx),
        lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE
      ),
      family   = binegbin_partialobs(),
      stanvars = binegbin_partialobs_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 20260706,
      refresh  = 0,
      init     = 0.5,
      control  = list(adapt_delta = 0.95)
    )
  })

  # Dispatch: both must resolve the package's log_lik_/posterior_predict_
  # methods by name without "no applicable method" errors.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))

  loo_obj <- suppressWarnings(brms::loo(fit))
  expect_s3_class(loo_obj, "loo")

  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))

  draws <- as.data.frame(fit)
  # Smoke gate at a wide interval, not a calibration claim. See
  # helper-coverage.R; the calibration assessment is in test-binegbin.R and
  # applies to this family's matched branch, which is binegbin term for term.
  check_recovery <- function(true_val, draws_col) {
    recovery_ok(draws, true_val, draws_col)
  }
  expect_true(check_recovery(true_log_mu_int,  "b_Intercept"))
  expect_true(check_recovery(log(true_lone),    "b_lamx_Intercept"))
  expect_true(check_recovery(true_sd_vessel,   "sd_vessel__Intercept"))
  expect_true(check_recovery(log(true_shapes), "b_shapes_Intercept"))
  expect_true(check_recovery(log(true_shapex), "b_shapexx_Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
