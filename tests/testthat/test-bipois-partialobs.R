# tests/testthat/test-bipois-partialobs.R
#
# The bipois likelihood on paired and unpaired rows: normalisation, the closed
# form of the unpaired branch, the Poisson limit of binegbin, conditional
# draws of y1 given y2, agreement of Stan with the R reference, and a brms fit
# of bipois_partialobs().

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("paired (y1_obs==1) joint PMF normalises to 1 across parameter sets", {
  norm_check <- function(mu, lone, ltwo, K = 100) {
    ys <- 0:K
    yg <- expand.grid(y1 = ys, y2 = ys)
    lp <- bipois_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo)
    sum(exp(lp))
  }
  params <- list(
    c(mu = 5,   lone = 3,  ltwo = 3),
    c(mu = 10,  lone = 2,  ltwo = 6),
    c(mu = 0.5, lone = 10, ltwo = 10)
  )
  for (p in params) {
    expect_equal(norm_check(p[["mu"]], p[["lone"]], p[["ltwo"]]),
                 1, tolerance = 1e-8, label = paste(p, collapse = ","))
  }
})

test_that("unpaired (y1_obs==0) branch normalises to 1 over y2", {
  # A Poisson pmf sums to 1 at any rate, so the sum checks that y1_obs == 0
  # selects a pmf in y2 alone. The assertions at the end check that the branch
  # does not depend on y1 or lambdaone.
  norm_check <- function(mu, lone, ltwo, K = 250) {
    ys <- 0:K
    lp <- bipois_lpmf_r(rep(0L, length(ys)), ys, y1_obs = 0L, mu, lone, ltwo)
    sum(exp(lp))
  }
  params <- list(
    c(mu = 5,  lone = 3, ltwo = 3),
    c(mu = 10, lone = 2, ltwo = 6),
    c(mu = 2,  lone = 4, ltwo = 4)
  )
  for (p in params) {
    expect_equal(norm_check(p[["mu"]], p[["lone"]], p[["ltwo"]]),
                 1, tolerance = 1e-8, label = paste(p, collapse = ","))
  }

  # The marginal of y2 does not depend on y1 or lambdaone.
  base <- bipois_lpmf_r(0L, 7L, y1_obs = 0L, mu = 5, lambdaone = 3, lambdatwo = 4)
  expect_equal(bipois_lpmf_r(99L, 7L, y1_obs = 0L, 5, 3, 4),  base)
  expect_equal(bipois_lpmf_r(0L,  7L, y1_obs = 0L, 5, 40, 4), base)
})

test_that("the closed form of the unpaired branch equals the brute-force convolution", {
  # With Poisson components, the convolution P(y2) = sum_k f_s(k) f_2(y2 - k)
  # equals Poisson(mu + lambdatwo), and the Stan and R code both use that
  # closed form. The brute-force sum below is the independent reference.
  convolve_r <- function(y2, mu, ltwo) {
    k <- 0:y2
    log_terms <- stats::dpois(k, mu, log = TRUE) +
      stats::dpois(y2 - k, ltwo, log = TRUE)
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }
  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    ltwo = c(0.2, 1, 5, 20)
  )
  y2_vals <- c(0L, 1L, 3L, 8L, 25L, 60L)
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    brute  <- vapply(y2_vals, convolve_r, numeric(1), mu = g$mu, ltwo = g$ltwo)
    closed <- bipois_lpmf_r(0L, y2_vals, y1_obs = 0L, g$mu, 1, g$ltwo)
    expect_equal(brute, closed, tolerance = 1e-12,
                 label = paste(unlist(g), collapse = ","))
  }
})

test_that("marginal identity: sum over y1 of the paired branch == unpaired branch", {
  # The paired branch (y1_obs == 1) summed over y1 must equal the unpaired
  # branch (y1_obs == 0) at the same y2.
  mu <- 6; lone <- 3; ltwo <- 4
  K  <- 200
  for (yl in c(0L, 1L, 3L, 7L, 15L)) {
    ys <- 0:K
    lp_joint <- bipois_lpmf_r(ys, rep(yl, length(ys)), mu, lone, ltwo)
    marg_from_joint <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
    marginal_branch <- bipois_lpmf_r(0L, yl, y1_obs = 0L, mu, lone, ltwo)
    expect_equal(marg_from_joint, marginal_branch, tolerance = 1e-10,
                 label = paste("y2 =", yl))
  }
})

# A call without y1_obs must evaluate the paired branch. A wrong default
# would score every row of a fully paired design against the marginal of y2,
# without an error.

test_that("a call without the indicator evaluates the paired branch, not the marginal", {
  grid <- expand.grid(
    mu   = c(0.5, 3, 12),
    lone = c(0.5, 2, 6),
    ltwo = c(0.5, 2, 6)
  )
  ys <- 0:20
  yg <- expand.grid(y1 = ys, y2 = ys)
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    default  <- bipois_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo)
    explicit <- bipois_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, y1_obs = 1L)
    expect_equal(default, explicit, tolerance = 1e-14,
                 label = paste(unlist(g), collapse = ","))
  }
  # The unpaired branch differs from the paired branch on more than 90% of
  # the grid, so the check above is not vacuous.
  g <- grid[1, ]
  marg <- bipois_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo, y1_obs = 0L)
  base <- bipois_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo)
  expect_gt(mean(abs(marg - base) > 1e-8), 0.9)
})

test_that("binegbin reduces to bipois in the Poisson limit, on both branches", {
  # The marginal identity for binegbin compares two sums from the same code.
  # As the three dispersions tend to infinity, binegbin tends to bipois, whose
  # unpaired branch is closed form, which gives an independent reference. The
  # error is of order 1/phi: 2.53e-3 at phi = 1e5 and 2.53e-5 at phi = 1e7 on
  # this grid. Each bound lies about 4 times above the observed error.
  mu <- 4; lone <- 2.5; ltwo <- 3
  ys <- 0:25
  yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))

  ref <- bipois_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, mu, lone, ltwo)

  err <- function(phi) {
    nb <- binegbin_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, mu, lone, ltwo,
                                shapes = phi, shapexone = phi, shapextwo = phi)
    max(abs(nb - ref))
  }

  e5 <- err(1e5)
  e7 <- err(1e7)
  expect_lt(e5, 1e-2)
  expect_lt(e7, 1e-4)
  # The error must decrease as phi increases.
  expect_lt(e7, e5)
})

test_that("conditional draws reproduce the joint/marginal conditional y1 | y2", {
  # .impute_draw_bipois() draws N_shared given y2 from a Binomial and adds a
  # new draw of the exclusive count N_1. The distribution of the draws must
  # equal P(y1 | y2) = joint(y1, y2) / marginal(y2), computed exactly. The
  # tolerances allow for Monte Carlo error with 2e5 draws.
  set.seed(20260804)
  mu <- 5; lone <- 3; ltwo <- 4
  y2 <- 6L
  ndraws <- 2e5

  prep <- make_synthetic_prep(
    dpars = list(
      mu        = rep(mu,   ndraws),
      lambdaone = rep(lone, ndraws),
      lambdatwo = rep(ltwo, ndraws)
    ),
    Y     = 0L,          # response value is unused by the conditional draw
    vint1 = y2,
    vint2 = 1L
  )
  draws <- .impute_draw_bipois(1, prep)
  expect_length(draws, ndraws)

  K  <- 80
  xs <- 0:K
  lp_joint <- bipois_lpmf_r(xs, rep(y2, length(xs)), mu, lone, ltwo)
  lp_marg  <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
  p_cond   <- exp(lp_joint - lp_marg)

  emp <- tabulate(draws + 1L, nbins = K + 1L) / ndraws
  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)

  # .impute_mean_bipois() computes E[y1 | y2] in closed form, so the
  # tolerance is 1e-10.
  ep <- .impute_mean_bipois(prep)
  expect_equal(dim(ep), c(ndraws, 1L))
  expect_equal(unique(as.vector(ep)), sum(xs * p_cond), tolerance = 1e-10)
})

test_that("the unpaired branch ignores y1_obs only where it should", {
  # The likelihood depends on y1_obs: on an unpaired row, the likelihood term
  # is the marginal probability of y2. The conditional expectation does not
  # depend on y1_obs, and is computed on every row.
  prep <- make_synthetic_prep(
    dpars = list(mu = 5, lambdaone = 3, lambdatwo = 4),
    Y     = 9L,
    vint1 = 6L,
    vint2 = 0L        # unpaired row
  )
  # On the unpaired row, log_lik is the Poisson(mu + lambdatwo) log density
  # of y2.
  expect_equal(
    log_lik_bipois(1, prep),
    stats::dpois(6L, 5 + 4, log = TRUE)
  )
  # The conditional expectation is the same on unpaired and paired rows.
  prep_paired <- prep
  prep_paired$data$vint2 <- 1L
  expect_equal(.impute_mean_bipois(prep),
               .impute_mean_bipois(prep_paired))
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

# bipois() and bipois_partialobs() call the same Stan lpmf, which is compared
# with the R reference on both branches.
stan_code <- paste0("functions {\n", bipois_stan_funs, "}\nmodel {}\n")

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

test_that("Stan bipois_lpmf matches the R reference (both branches)", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5)
  )

  check_one <- function(mu, lone, ltwo) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
    stan_vals <- mapply(
      function(r, s, e) bipois_lpmf(r, mu, lone, ltwo, s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- bipois_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs, mu, lone, ltwo)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("the indicator selects the same branch in Stan and in the R reference", {
  # The observation indicator is the last argument of the Stan lpmf and a
  # named argument of the R reference. The branch selected by each value of
  # the indicator must agree between Stan and R.
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5)
  )
  check_one <- function(mu, lone, ltwo, flag) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(function(r, s) bipois_lpmf(r, mu, lone, ltwo, s, flag),
                        yg$y1, yg$y2)
    r_vals <- bipois_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, y1_obs = flag)
    max(abs(stan_vals - r_vals))
  }
  for (flag in c(1L, 0L)) {
    diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, MoreArgs = list(flag = flag))
    expect_true(max(diffs) < 1e-12,
                label = paste("y1_obs =", flag, "max diff =", max(diffs)))
  }

  # The two branches differ by more than 0.5 on the log scale, so an indicator
  # that selected the wrong branch would fail the loop above.
  paired <- bipois_lpmf(5L, 5, 3, 4, 6L, 1L)
  marg    <- bipois_lpmf(5L, 5, 3, 4, 6L, 0L)
  expect_gt(abs(paired - marg), 0.5)
})

test_that("Stan bipois_lpmf is numerically stable at extreme rates", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  edge_cases <- list(
    c(mu = 1e-6, lone = 1e-6, ltwo = 1e-6),
    c(mu = 1e-6, lone = 50,   ltwo = 50),
    c(mu = 200,  lone = 1e-6, ltwo = 1e-6),
    c(mu = 200,  lone = 200,  ltwo = 200)
  )
  for (ec in edge_cases) {
    ys <- 0:5
    yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))
    stan_vals <- mapply(
      function(r, s, e) bipois_lpmf(r, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- bipois_lpmf_r(yg$y1, yg$y2, y1_obs = yg$y1_obs,
                                  ec[["mu"]], ec[["lone"]], ec[["ltwo"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# brms end-to-end: dispatch (loo + posterior_predict) and recovery
# -----------------------------------------------------------------------

test_that("bipois_partialobs fits a partially paired design, dispatches, and recovers params", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260804)

  n_group     <- 8L
  n_per_group <- 30L
  n            <- n_group * n_per_group

  true_log_mu_int <- log(8)
  true_sd_group  <- 0.3
  true_lone       <- 3           # lambdaone = lambdatwo

  group <- rep(seq_len(n_group), each = n_per_group)
  z_mu   <- rnorm(n_group)
  mu_i   <- exp(true_log_mu_int + true_sd_group * z_mu[group])

  n_shared <- rpois(n, mu_i)
  n1       <- rpois(n, true_lone)
  n2       <- rpois(n, true_lone)

  # Half the rows are unpaired. The likelihood term for an unpaired row
  # depends on mu and lambdatwo only through mu + lambdatwo, and does not
  # depend on lambdaone, so only the paired rows separate mu, lambdaone, and
  # lambdatwo; see ?bipois_partialobs.
  y1_obs <- rep(c(1L, 0L), length.out = n)

  dat <- data.frame(
    y1     = n_shared + n1,
    y2     = n_shared + n2,
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
        lamx ~ 1, nl = TRUE
      ),
      family   = bipois_partialobs(),
      stanvars = bipois_partialobs_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 20260804,
      refresh  = 0,
      init     = 0.5,
      control  = list(adapt_delta = 0.95)
    )
  })

  # brms must find the log_lik_bipois() and posterior_predict_bipois() methods
  # of this package.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))

  loo_obj <- suppressWarnings(brms::loo(fit))
  expect_s3_class(loo_obj, "loo")

  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))

  ep <- posterior_epred_bipois(brms::prepare_predictions(fit))
  expect_equal(dim(ep)[2], n)
  expect_true(all(is.finite(ep)))

  # For an untruncated fit, brms::posterior_epred() dispatches to the family
  # method, so the two calls must agree exactly.
  expect_equal(brms::posterior_epred(fit), ep)

  draws <- as.data.frame(fit)
  # Smoke gates at a wide interval; see helper-coverage.R.
  expect_true(recovery_ok(draws, true_log_mu_int, "b_Intercept"))
  expect_true(recovery_ok(draws, log(true_lone),  "b_lamx_Intercept"))
  expect_true(recovery_ok(draws, true_sd_group,  "sd_group__Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
