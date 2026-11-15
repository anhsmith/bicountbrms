# tests/testthat/test-binegbin-partialobs.R
#
# The binegbin likelihood on paired and unpaired rows: normalisation, the
# marginal identity, the default branch, conditional draws of y1 given y2,
# agreement of the Stan lpmf with the R reference, and a brms fit of
# binegbin_partialobs().

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("paired (y1_obs==1) joint PMF normalises to 1 across parameter sets", {
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

test_that("the unpaired branch (y1_obs == 0) normalises to 1 over y2", {
  # The unpaired branch (y1_obs == 0) is a pmf in y2 alone, because the y1
  # margin is summed out.
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

test_that("marginal identity: the sum over y1 of the paired branch equals the unpaired branch", {
  # The paired branch (y1_obs == 1) summed over y1 must equal the unpaired
  # branch (y1_obs == 0) at the same y2.
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

# A call without y1_obs must evaluate the paired branch. A wrong default
# would score every row of a fully paired design against the marginal of y2.

test_that("a call without the indicator evaluates the paired branch, not the marginal", {
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
  # The unpaired branch differs from the paired branch on more than 90% of
  # the grid, so the check above is not vacuous.
  g <- grid[1, ]
  marg <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = 0L,
                          g$mu, g$lone, g$ltwo, g$ss, g$sx)
  base <- binegbin_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, g$ss, g$sx)
  expect_gt(mean(abs(marg - base) > 1e-8), 0.9)
})

test_that("conditional draws reproduce the joint/marginal conditional y1 | y2", {
  # .impute_draw_binegbin() draws N_shared given y2 and adds a new draw of the
  # exclusive count N_1. The distribution of the draws must equal
  # P(y1 | y2) = joint(y1, y2) / marginal(y2), computed exactly. The
  # tolerances allow for Monte Carlo error with 2e5 draws.
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
    Y     = 0L,          # response value is unused by the conditional draw
    vint1 = y2,
    vint2 = 1L
  )
  draws <- .impute_draw_binegbin(1, prep)
  expect_length(draws, ndraws)

  # Analytic conditional P(y1 = x | y2) = joint(x, y2) / marginal(y2).
  K <- 80
  xs <- 0:K
  lp_joint <- binegbin_lpmf_r(xs, rep(y2, length(xs)), mu, lone, ltwo, ss, sx)
  lp_marg  <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
  p_cond   <- exp(lp_joint - lp_marg)

  emp <- tabulate(draws + 1L, nbins = K + 1L) / ndraws
  # The comparison is restricted to values with probability above 1e-3.
  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

# binegbin() and binegbin_partialobs() call the same Stan lpmf, which is
# compared with the R reference on both branches.
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
    # The single `sx` is passed for both shapexone and shapextwo. The grid
    # with unequal dispersions is in test-binegbin-dispersions.R.
    stan_vals <- mapply(
      function(r, s, e) binegbin_lpmf(r, mu, lone, ltwo, ss, sx, sx, s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, mu, lone, ltwo, ss, sx)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, grid$ss, grid$sx)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("the Stan paired branch routes the two dispersions the way R does", {
  # The previous test passes one value for both dispersions and so cannot
  # detect an exchange of shapexone and shapextwo between the Stan signature
  # and the R reference. This test sets shapexone = 0.6 and shapextwo = 7.
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5),
    ss   = c(0.5, 2, 50)
  )
  SX1 <- 0.6; SX2 <- 7      # an order of magnitude apart

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

  # Exchanging the two dispersions in the Stan call changes the log density by
  # more than 0.05, so the agreement above constrains the argument order.
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
  # The simulation uses one dispersion for both exclusive components. The
  # model ties shapexone and shapextwo to the single non-linear parameter
  # `shapexx`.
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260706)

  n_group     <- 8L
  n_per_group <- 30L
  n            <- n_group * n_per_group

  true_log_mu_int <- log(8)
  true_sd_group  <- 0.3
  true_lone        <- 3           # lambdaone = lambdatwo
  true_shapes     <- 2
  true_shapex     <- 1.5

  group <- rep(seq_len(n_group), each = n_per_group)
  z_mu   <- rnorm(n_group)
  mu_i   <- exp(true_log_mu_int + true_sd_group * z_mu[group])

  n_shared <- rnbinom(n, size = true_shapes, mu = mu_i)
  n1      <- rnbinom(n, size = true_shapex, mu = true_lone)
  n2      <- rnbinom(n, size = true_shapex, mu = true_lone)

  # Half the rows are unpaired.
  y1_obs <- rep(c(1L, 0L), length.out = n)

  dat <- data.frame(
    y1   = n_shared + n1,
    y2   = n_shared + n2,
    y1_obs = y1_obs,
    group = factor(group)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2, y1_obs) ~ 1,
        mu ~ 1 + (1 | group),
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

  # brms must find the log_lik_binegbin() and posterior_predict_binegbin()
  # methods of this package.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))

  loo_obj <- suppressWarnings(brms::loo(fit))
  expect_s3_class(loo_obj, "loo")

  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))

  draws <- as.data.frame(fit)
  # Smoke gates at a wide interval; see helper-coverage.R. The calibration
  # assessment of the paired branch is in test-binegbin.R.
  check_recovery <- function(true_val, draws_col) {
    recovery_ok(draws, true_val, draws_col)
  }
  expect_true(check_recovery(true_log_mu_int,  "b_Intercept"))
  expect_true(check_recovery(log(true_lone),    "b_lamx_Intercept"))
  expect_true(check_recovery(true_sd_group,   "sd_group__Intercept"))
  expect_true(check_recovery(log(true_shapes), "b_shapes_Intercept"))
  expect_true(check_recovery(log(true_shapex), "b_shapexx_Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
