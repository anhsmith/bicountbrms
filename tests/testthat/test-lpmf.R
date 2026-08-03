# tests/testthat/test-lpmf.R

# R-side log-PMF: same computation as log_lik_skellam1, now expressed in
# terms of sigma (the sampled parameter); mu = sigma^2 / 2 is the
# underlying Skellam(mu, mu) rate fed to the unchanged Bessel formula.
r_lpmf <- function(k, sigma) {
  mu <- sigma^2 / 2
  log(besselI(2 * mu, abs(k), expon.scaled = TRUE))
}

# mu grid carried over unchanged from the pre-reparameterisation tests;
# sigma = sqrt(2 * mu) is the corresponding sampled-parameter grid.
mu_vals    <- c(0.5, 1, 2, 5, 10, 20)
sigma_vals <- sqrt(2 * mu_vals)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF matches dskellam across grid", {
  k_vals <- -10L:10L

  for (i in seq_along(mu_vals)) {
    mu    <- mu_vals[i]
    sigma <- sigma_vals[i]
    r_vals   <- vapply(k_vals, r_lpmf, numeric(1), sigma = sigma)
    ref_vals <- skellam::dskellam(k_vals, lambda1 = mu, lambda2 = mu, log = TRUE)
    expect_equal(r_vals, ref_vals, tolerance = 1e-6,
                 label = paste0("mu = ", mu))
  }
})

test_that("R log-PMF handles k = 0 correctly", {
  r_vals   <- vapply(sigma_vals, r_lpmf, numeric(1), k = 0L)
  ref_vals <- skellam::dskellam(0L, lambda1 = mu_vals, lambda2 = mu_vals, log = TRUE)
  expect_equal(r_vals, ref_vals, tolerance = 1e-6)
})

test_that("R log-PMF is numerically stable at large sigma", {
  for (mu in c(100, 500)) {
    sigma <- sqrt(2 * mu)
    for (k in c(-5L, 0L, 5L)) {
      val <- r_lpmf(k, sigma)
      expect_false(is.nan(val),      label = paste0("NaN at sigma=", sigma, " k=", k))
      expect_false(is.infinite(val), label = paste0("Inf at sigma=", sigma, " k=", k))
    }
  }
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

stan_code <- paste0("functions {\n", skellam1_stan_funs, "}\nmodel {}\n")

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

test_that("Stan log-PMF matches dskellam across grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  k_vals <- -10L:10L

  for (i in seq_along(mu_vals)) {
    mu    <- mu_vals[i]
    sigma <- sigma_vals[i]
    stan_vals <- vapply(k_vals, function(k) skellam1_lpmf(k, sigma), numeric(1))
    ref_vals  <- skellam::dskellam(k_vals, lambda1 = mu, lambda2 = mu, log = TRUE)
    expect_equal(stan_vals, ref_vals, tolerance = 1e-6,
                 label = paste0("mu = ", mu))
  }
})

test_that("Stan log-PMF matches R log-PMF across grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  k_vals <- -10L:10L

  for (sigma in sigma_vals) {
    stan_vals <- vapply(k_vals, function(k) skellam1_lpmf(k, sigma), numeric(1))
    r_vals    <- vapply(k_vals, r_lpmf, numeric(1), sigma = sigma)
    expect_equal(stan_vals, r_vals, tolerance = 1e-6,
                 label = paste0("sigma = ", sigma))
  }
})

test_that("Stan log-PMF is numerically stable at large sigma", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  for (mu in c(100, 500)) {
    sigma <- sqrt(2 * mu)
    for (k in c(-5L, 0L, 5L)) {
      val <- skellam1_lpmf(k, sigma)
      expect_false(is.nan(val),      label = paste0("NaN at sigma=", sigma, " k=", k))
      expect_false(is.infinite(val), label = paste0("Inf at sigma=", sigma, " k=", k))
    }
  }
})
