# tests/testthat/test-bipois.R

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force lpmf normalises to 1 across parameter sets", {
  norm_check <- function(mu, lambdaone, lambdatwo, K = 60) {
    ys <- 0:K
    yg <- expand.grid(y1 = ys, y2 = ys)
    lp <- bipois_lpmf_r(yg$y1, yg$y2, mu, lambdaone, lambdatwo)
    sum(exp(lp))
  }
  for (p in list(c(1, 1, 1), c(5, 2, 3), c(0.5, 10, 10))) {
    expect_equal(norm_check(p[1], p[2], p[3]), 1, tolerance = 1e-8,
                 label = paste(p, collapse = ","))
  }
})

test_that("R brute-force lpmf reproduces Poisson(mu+lambdaone) marginal for y1", {
  mu <- 3; lambdaone <- 2; lambdatwo <- 4
  K <- 80
  marg_em <- vapply(0:K, function(r) {
    s_vals <- 0:K
    lp <- bipois_lpmf_r(rep(r, length(s_vals)), s_vals, mu, lambdaone, lambdatwo)
    mx <- max(lp)
    mx + log(sum(exp(lp - mx)))
  }, numeric(1))
  ref_em <- dpois(0:K, mu + lambdaone, log = TRUE)
  expect_equal(marg_em, ref_em, tolerance = 1e-6)
})

test_that("N_shared cancels: d = y1 - y2 marginal matches Skellam(lambdaone, lambdatwo)", {
  # `skellam` is a Suggests, not an Imports, since 0.9.0 -- no function in R/
  # calls it. It is still the right independent reference for this identity,
  # so the test guards rather than dropping it.
  skip_if_not_installed("skellam")
  mu <- 4; lambdaone <- 2; lambdatwo <- 3
  K <- 60
  ys <- 0:K
  yg <- expand.grid(y1 = ys, y2 = ys)
  lp <- bipois_lpmf_r(yg$y1, yg$y2, mu, lambdaone, lambdatwo)
  d  <- yg$y1 - yg$y2
  d_vals <- -10:10
  p_d <- vapply(d_vals, function(dd) sum(exp(lp[d == dd])), numeric(1))
  ref <- skellam::dskellam(d_vals, lambda1 = lambdaone, lambda2 = lambdatwo)
  expect_equal(p_d, ref, tolerance = 1e-4)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

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

test_that("Stan bipois_lpmf matches R brute-force reference across a grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu       = c(1e-6, 0.1, 1, 5, 20, 100),
    lambdaone = c(1e-6, 0.1, 1, 5, 20, 100),
    lambdatwo = c(1e-6, 0.1, 1, 5, 20, 100)
  )

  check_one <- function(mu, lambdaone, lambdatwo) {
    means <- c(mu + lambdaone, mu + lambdatwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 5L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(function(r, s) bipois_lpmf(r, mu, lambdaone, lambdatwo, s, 1L),
                         yg$y1, yg$y2)
    r_vals <- bipois_lpmf_r(yg$y1, yg$y2, mu, lambdaone, lambdatwo)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lambdaone, grid$lambdatwo)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("Stan bipois_lpmf is numerically stable at near-zero and large rates", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  edge_cases <- list(
    c(mu = 1e-8, lambdaone = 1e-8, lambdatwo = 1e-8),
    c(mu = 1e-8, lambdaone = 50,   lambdatwo = 50),
    c(mu = 200,  lambdaone = 1e-8, lambdatwo = 1e-8),
    c(mu = 200,  lambdaone = 200,  lambdatwo = 200)
  )

  for (ec in edge_cases) {
    ys <- 0:5
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(function(r, s) bipois_lpmf(r, ec[["mu"]], ec[["lambdaone"]], ec[["lambdatwo"]], s, 1L),
                         yg$y1, yg$y2)
    r_vals <- bipois_lpmf_r(yg$y1, yg$y2, ec[["mu"]], ec[["lambdaone"]], ec[["lambdatwo"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# Parameter recovery from simulated hierarchical data (brms end-to-end)
# -----------------------------------------------------------------------

test_that("bipois parameter recovery from simulated vessel-level data", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260704)

  n_vessel     <- 8L
  n_per_vessel <- 25L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int       <- log(3)
  true_log_lambdaone_int <- log(1.5)
  true_log_lambdatwo_int <- log(4)
  true_sd_vessel        <- 0.4

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)

  z_mu       <- rnorm(n_vessel)
  z_lambdaone <- rnorm(n_vessel)
  z_lambdatwo <- rnorm(n_vessel)

  mu_i       <- exp(true_log_mu_int       + true_sd_vessel * z_mu[vessel])
  lambdaone_i <- exp(true_log_lambdaone_int + true_sd_vessel * z_lambdaone[vessel])
  lambdatwo_i <- exp(true_log_lambdatwo_int + true_sd_vessel * z_lambdatwo[vessel])

  n_shared <- rpois(n, mu_i)
  n1      <- rpois(n, lambdaone_i)
  n2      <- rpois(n, lambdatwo_i)

  dat <- data.frame(
    y1   = n_shared + n1,
    y2   = n_shared + n2,
    vessel = factor(vessel)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2) ~ 1,
        mu       ~ 1 + (1 | vessel),
        lambdaone ~ 1 + (1 | vessel),
        lambdatwo ~ 1 + (1 | vessel)
      ),
      family   = bipois(),
      stanvars = bipois_stanvars(),
      data     = dat,
      backend  = "rstan",
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 20260704,
      refresh  = 0,
      control  = list(adapt_delta = 0.95)
    )
  })

  draws <- as.data.frame(fit)

  # Smoke gate at a wide interval, not a calibration claim: a 90% interval
  # fails 10% of the time on a CORRECT model, which across the six assertions
  # below is a ~46% chance of a spurious failure. See helper-coverage.R.
  check_recovery <- function(true_val, draws_col) {
    recovery_ok(draws, true_val, draws_col)
  }

  # brms treats "mu" as the family's canonical/default dpar and drops its
  # infix from generated column names (b_Intercept, sd_vessel__Intercept),
  # unlike the other two, plainly-named dpars (b_lambdaone_Intercept, etc.)
  # -- so "mu" here reads as b_Intercept even though it stands in for the
  # shared rate rather than any mean.
  expect_true(check_recovery(true_log_mu_int,       "b_Intercept"))
  expect_true(check_recovery(true_log_lambdaone_int, "b_lambdaone_Intercept"))
  expect_true(check_recovery(true_log_lambdatwo_int, "b_lambdatwo_Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__lambdaone_Intercept"))
  expect_true(check_recovery(true_sd_vessel, "sd_vessel__lambdatwo_Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
