# tests/testthat/test-dnorm2.R
#
# Validates dnorm2(), the free-location/free-scale discrete normal
# family, with the same rigor as dnorm1's test suite. See
# test-dnorm1.R's file header for the full citation of the documented
# Stan limitation (stan-dev/math#1985: normal_lccdf underflows to -inf
# above (y-mu)/sigma ~ 8.25) and why the erfc()-based exact survival
# form is used in stanfunctions.R instead -- both apply identically here.

log_normal_cdf   <- function(x, sigma) pnorm(x, sd = sigma, log.p = TRUE)
log_normal_lccdf <- function(x, sigma) pnorm(x, sd = sigma, lower.tail = FALSE, log.p = TRUE)
log_diff_exp_r   <- function(a, bb) a + log1p(-exp(bb - a))

# R-side log-PMF: same stable log-space construction as test-dnorm1.R's
# reference, branched on whether z is on the far side of mu rather than
# of 0.
r_lpmf_dnorm2 <- function(z, mu, sigma) {
  ifelse(
    z >= mu,
    log_diff_exp_r(log_normal_lccdf(z - mu - 0.5, sigma), log_normal_lccdf(z - mu + 0.5, sigma)),
    log_diff_exp_r(log_normal_cdf(z - mu + 0.5, sigma), log_normal_cdf(z - mu - 0.5, sigma))
  )
}

mu_vals    <- c(-20, -5, 0, 5, 20)
sigma_vals <- c(0.5, 1, 5, 20, 50)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF normalises to 1 across grid, including nonzero mu", {
  for (mu in mu_vals) {
    for (sigma in sigma_vals) {
      # kmax anchored to ~10 SDs, not dlaplace2's 40*sigma+50 -- see
      # test-dnorm1.R for why the discrete normal needs far fewer SDs to
      # capture essentially all the mass than the heavier-tailed Laplace.
      half  <- ceiling(10 * sigma + 20)
      total <- sum(exp(r_lpmf_dnorm2((mu - half):(mu + half), mu, sigma)))
      expect_equal(total, 1, tolerance = 1e-8,
                   label = paste0("mu = ", mu, ", sigma = ", sigma))
    }
  }
})

test_that("mu = 0 reduces exactly to dnorm1's lpmf", {
  r_lpmf_dnorm1 <- function(z, sigma) {
    ifelse(
      z >= 0,
      log_diff_exp_r(log_normal_lccdf(z - 0.5, sigma), log_normal_lccdf(z + 0.5, sigma)),
      log_diff_exp_r(log_normal_cdf(z + 0.5, sigma), log_normal_cdf(z - 0.5, sigma))
    )
  }
  for (sigma in c(0.5, 5, 20)) {
    for (z in -5:5) {
      expect_equal(r_lpmf_dnorm2(z, 0, sigma), r_lpmf_dnorm1(z, sigma), tolerance = 1e-10,
                   label = paste0("sigma = ", sigma, ", z = ", z))
    }
  }
})

test_that("mu shifts the distribution exactly: r_lpmf_dnorm2(z, mu, sigma) == r_lpmf_dnorm2(z - mu, 0, sigma)", {
  for (mu in mu_vals) {
    for (sigma in c(1, 10)) {
      for (z_offset in -5:5) {
        z <- round(mu) + z_offset
        expect_equal(r_lpmf_dnorm2(z, mu, sigma), r_lpmf_dnorm2(z - mu, 0, sigma), tolerance = 1e-10,
                     label = paste0("mu = ", mu, ", sigma = ", sigma, ", z = ", z))
      }
    }
  }
})

test_that("mu and sigma are not coupled: small sigma with large |mu| is valid (unlike skellam2)", {
  mu <- 1000
  sigma <- 0.5
  total <- sum(exp(r_lpmf_dnorm2((mu - 50):(mu + 50), mu, sigma)))
  expect_equal(total, 1, tolerance = 1e-6)
})

test_that("R log-PMF is numerically stable, realistic-but-stressed range", {
  for (mu in c(-50, 0, 50)) {
    for (sigma in c(0.2, 1, 10, 100)) {
      for (k in round(mu + sigma * c(-10, 0, 10))) {
        val <- r_lpmf_dnorm2(k, mu, sigma)
        expect_false(is.nan(val),      label = paste0("NaN at mu=", mu, " sigma=", sigma, " k=", k))
        expect_false(is.infinite(val), label = paste0("Inf at mu=", mu, " sigma=", sigma, " k=", k))
      }
    }
  }
})

log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Brute-force log P(Z > y), summing from y+1 (dnorm2_lccdf is defined
# directly as log P(Z>y), no skellam-style offset -- see test-dnorm1.R).
# K anchored to max(y, mu), not y alone: when y is far *below* mu (e.g.
# probing -10 SDs from the mean), the sum still needs to extend ~10 SDs
# *above* mu to capture the upper tail that makes up most of P(Z>y) in
# that regime -- anchoring only to y undershoots badly there (confirmed:
# anchoring to y alone produced a ~0.42 log-probability error at
# mu=20, sigma=50, y_mult=-10, where K landed at only mu+20, nowhere
# near far enough past the mean for sigma=50).
logS_bruteforce_dnorm2 <- function(y, mu, sigma, K = NULL) {
  if (is.null(K)) K <- max(y, mu) + ceiling(10 * sigma + 20)
  log_sum_exp(r_lpmf_dnorm2((y + 1):K, mu, sigma))
}

# Closed-form log P(Z > y) via the same stable survival function used in
# the reference's z>=mu branch.
logS_closed_dnorm2 <- function(y, mu, sigma) log_normal_lccdf(y - mu + 0.5, sigma)

y_multiples <- c(-10, -5, -3, -1, 0, 1, 3, 5, 10)
trunc_grid  <- expand.grid(mu = mu_vals, sigma = sigma_vals, y_mult = y_multiples)
trunc_grid$y <- round(trunc_grid$mu + trunc_grid$sigma * trunc_grid$y_mult)

test_that("R brute-force log-CCDF agrees with the closed form across grid", {
  diffs <- mapply(
    function(mu, sigma, y) logS_closed_dnorm2(y, mu, sigma) - logS_bruteforce_dnorm2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-8)
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lpmf_stan_code <- paste0("functions {\n", dnorm2_stan_funs, "}\nmodel {}\n")

stan_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lpmf_stan_code)
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan log-PMF matches R log-space reference across grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  for (mu in mu_vals) {
    for (sigma in sigma_vals) {
      k_vals    <- unique(round(mu + sigma * (-10:10)))
      stan_vals <- vapply(k_vals, function(k) dnorm2_lpmf(as.integer(k), mu, sigma), numeric(1))
      r_vals    <- vapply(k_vals, r_lpmf_dnorm2, numeric(1), mu = mu, sigma = sigma)
      expect_equal(stan_vals, r_vals, tolerance = 1e-6,
                   label = paste0("mu = ", mu, ", sigma = ", sigma))
    }
  }
})

lccdf_stan_code <- paste0("functions {\n", dnorm2_lccdf_stan, "}\nmodel {}\n")

lccdf_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm_lccdf <- rstan::stan_model(model_code = lccdf_stan_code)
      rstan::expose_stan_functions(sm_lccdf)
    })
    lccdf_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan dnorm2_lccdf matches R closed-form / brute-force references", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  diffs_closed <- mapply(
    function(mu, sigma, y) dnorm2_lccdf(as.integer(y), mu, sigma) - logS_closed_dnorm2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_closed)), 1e-6)

  diffs_brute <- mapply(
    function(mu, sigma, y) dnorm2_lccdf(as.integer(y), mu, sigma) - logS_bruteforce_dnorm2(y, mu, sigma),
    trunc_grid$mu, trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_brute)), 1e-6)
})

test_that("Stan dnorm2_lpmf/lccdf are numerically stable, realistic-but-stressed range", {
  skip_if_not(stan_ready && lccdf_ready, "rstan unavailable or Stan compilation failed")
  for (mu in c(-50, 0, 50)) {
    for (sigma in c(0.2, 1, 10, 100)) {
      for (k in round(mu + sigma * c(-10, 0, 10))) {
        lpmf_val  <- dnorm2_lpmf(as.integer(k), mu, sigma)
        lccdf_val <- dnorm2_lccdf(as.integer(k), mu, sigma)
        expect_false(is.nan(lpmf_val) || is.infinite(lpmf_val),
                     label = paste0("lpmf at mu=", mu, " sigma=", sigma, " k=", k))
        expect_false(is.nan(lccdf_val) || is.infinite(lccdf_val),
                     label = paste0("lccdf at mu=", mu, " sigma=", sigma, " k=", k))
      }
    }
  }
})

# -----------------------------------------------------------------------
# Structural check: no constraint coupling mu and sigma
# -----------------------------------------------------------------------

test_that("make_stancode() shows no constraint coupling mu and sigma", {
  skip_if_not_installed("brms")

  fake_dat <- data.frame(delta = sample(-5:5, 20, replace = TRUE), x = rnorm(20))
  code <- brms::make_stancode(
    brms::bf(delta ~ x, sigma ~ x),
    family   = dnorm2(),
    stanvars = dnorm2_stanvars(),
    data     = fake_dat
  )
  expect_false(grepl("reject\\(", code))
})

# -----------------------------------------------------------------------
# End-to-end smoke tests: resp_trunc() recovery, zero and nonzero mu_true
# -----------------------------------------------------------------------

fit_and_check_dnorm2 <- function(mu_true, sigma_true, seed) {
  set.seed(seed)
  n    <- 60
  y2 <- sample(0:10, n, replace = TRUE)

  draw_trunc <- function(y) {
    repeat {
      d <- round(rnorm(1, mean = mu_true, sd = sigma_true))
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y2, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y2)

  sane_prior <- brms::prior(normal(0, 10), class = "Intercept") +
    brms::prior(normal(1, 1.5), class = "Intercept", dpar = "sigma")

  fit <- brms::brm(
    brms::bf(delta | trunc(lb = neg_bound) ~ 1, sigma ~ 1),
    family   = dnorm2(),
    stanvars = dnorm2_stanvars() + dnorm2_lccdf_stanvars(),
    data     = dat,
    prior    = sane_prior,
    backend  = "cmdstanr",
    chains   = 2,
    iter     = 1000,
    warmup   = 500,
    seed     = seed,
    refresh  = 0
  )

  draws     <- as.data.frame(fit)
  mu_hat    <- draws[["b_Intercept"]]
  mu_q      <- quantile(mu_hat, c(0.025, 0.975))
  sigma_hat <- exp(draws[["b_sigma_Intercept"]])
  sigma_q   <- quantile(sigma_hat, c(0.025, 0.975))

  list(
    fit  = fit,
    dat  = dat,
    mu_ok = mu_true >= mu_q[[1]] && mu_true <= mu_q[[2]],
    mu_label = paste0("mu_true = ", mu_true, ", 95% CI: [", round(mu_q[[1]], 3), ", ", round(mu_q[[2]], 3), "]"),
    sigma_ok = sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    sigma_label = paste0("sigma_true = ", sigma_true, ", 95% CI: [", round(sigma_q[[1]], 3), ", ", round(sigma_q[[2]], 3), "]")
  )
}

test_that("dnorm2 recovers a genuinely nonzero mu_true under resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  suppressMessages({
    res <- fit_and_check_dnorm2(mu_true = 4, sigma_true = 5, seed = 8)
  })

  expect_s3_class(res$fit, "brmsfit")
  expect_true(res$mu_ok, label = res$mu_label)
  expect_true(res$sigma_ok, label = res$sigma_label)

  # Extend the existing fit (no new brm() call) to also check that
  # posterior_predict()/posterior_epred() honour resp_trunc()'s bounds --
  # confirms the R-side truncation fix works through the full brms
  # dispatch path, not just the synthetic-prep unit tests below.
  pp <- brms::posterior_predict(res$fit)
  expect_true(all(sweep(pp, 2, res$dat$neg_bound, `>=`)),
              label = "posterior_predict draws below their row's trunc(lb=) bound")

  # brms::posterior_epred() itself cannot be used here -- see the matching
  # comment in test-lccdf.R for the confirmed brms-level reason
  # (posterior_epred.brmsprep() unconditionally routes truncated fits to
  # posterior_epred_trunc(), which has no generic handler for ANY custom
  # family and always errors, independent of this package's own fix).
  prep <- brms::prepare_predictions(res$fit)
  ep   <- posterior_epred_dnorm2(prep)
  expect_true(all(is.finite(ep)))
  expect_true(all(sweep(ep, 2, res$dat$neg_bound, `>=`)),
              label = "posterior_epred below the row's trunc(lb=) bound")
})

test_that("dnorm2 does not spuriously detect bias when mu_true = 0", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  suppressMessages({
    res <- fit_and_check_dnorm2(mu_true = 0, sigma_true = 5, seed = 9)
  })

  expect_s3_class(res$fit, "brmsfit")
  expect_true(res$mu_ok, label = res$mu_label)
  expect_true(res$sigma_ok, label = res$sigma_label)

  pp <- brms::posterior_predict(res$fit)
  expect_true(all(sweep(pp, 2, res$dat$neg_bound, `>=`)),
              label = "posterior_predict draws below their row's trunc(lb=) bound")

  # brms::posterior_epred() itself cannot be used here -- see the matching
  # comment in test-lccdf.R for the confirmed brms-level reason
  # (posterior_epred.brmsprep() unconditionally routes truncated fits to
  # posterior_epred_trunc(), which has no generic handler for ANY custom
  # family and always errors, independent of this package's own fix).
  prep <- brms::prepare_predictions(res$fit)
  ep   <- posterior_epred_dnorm2(prep)
  expect_true(all(is.finite(ep)))
  expect_true(all(sweep(ep, 2, res$dat$neg_bound, `>=`)),
              label = "posterior_epred below the row's trunc(lb=) bound")
})

# -----------------------------------------------------------------------
# posterior_predict / posterior_epred truncation correctness (synthetic prep)
# -----------------------------------------------------------------------

test_that("posterior_predict_dnorm2 respects lb (repro of confirmed bug)", {
  # mu=1.5, sigma=4, lb=-14: the exact combination confirmed (during manual
  # testing) to produce out-of-bound draws under the pre-fix code -- range
  # [-20, 17] observed there, vs. the correct [-14, ...] here.
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(1.5, nrow = 3000, ncol = 1),
                 sigma = matrix(4, nrow = 3000, ncol = 1)),
    Y = 0, lb = -14
  )
  set.seed(123)
  draws <- posterior_predict_dnorm2(1, prep)
  expect_true(all(draws >= -14), label = paste0("min draw = ", min(draws)))
})

test_that("posterior_predict_dnorm2 without lb/ub matches untruncated behaviour (fast path)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(1.5, nrow = 5000, ncol = 1),
                 sigma = matrix(4, nrow = 5000, ncol = 1)),
    Y = 0
  )
  set.seed(1)
  draws <- posterior_predict_dnorm2(1, prep)
  expect_true(is.numeric(draws) && length(draws) == 5000)
  expect_true(min(draws) < -5 && max(draws) > 8)
})

test_that("posterior_predict_dnorm2 draws match the truncated PMF distributionally", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(1.5, nrow = 20000, ncol = 1),
                 sigma = matrix(4, nrow = 20000, ncol = 1)),
    Y = 0, lb = -5, ub = 8
  )
  set.seed(42)
  draws <- posterior_predict_dnorm2(1, prep)
  expect_true(all(draws >= -5 & draws <= 8))

  support <- -5:8
  probs <- exp(dnorm2_lpmf_r(support, 1.5, 4))
  probs <- probs / sum(probs)
  emp <- as.numeric(table(factor(draws, levels = support))) / length(draws)
  expect_lt(max(abs(emp - probs)), 0.02)
})

test_that("posterior_epred_dnorm2 differs from untruncated mu when lb is tight", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(1.5, nrow = 5, ncol = 1),
                 sigma = matrix(4, nrow = 5, ncol = 1)),
    Y = 0, lb = -2
  )
  epred <- posterior_epred_dnorm2(prep)
  expect_equal(dim(epred), c(5, 1))
  expect_true(all(epred[, 1] > 1.5), label = paste0("epred = ", epred[1, 1]))
})

test_that("posterior_epred_dnorm2 matches brute-force truncated mean", {
  mu_val <- 1.5; sigma_val <- 4; lb_val <- -14
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(mu_val, nrow = 1, ncol = 1),
                 sigma = matrix(sigma_val, nrow = 1, ncol = 1)),
    Y = 0, lb = lb_val
  )
  epred <- posterior_epred_dnorm2(prep)

  support <- lb_val:round(mu_val + 10 * sigma_val)
  probs <- exp(dnorm2_lpmf_r(support, mu_val, sigma_val))
  probs <- probs / sum(probs)
  brute_force_mean <- sum(support * probs)

  expect_equal(epred[1, 1], brute_force_mean, tolerance = 1e-6)
})

test_that("posterior_epred_dnorm2 leaves untruncated observations exactly at mu (no regression)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(c(1.5, -3), nrow = 2, ncol = 2),
                 sigma = matrix(4, nrow = 2, ncol = 2)),
    Y = c(0, 0)
  )
  epred <- posterior_epred_dnorm2(prep)
  expect_equal(epred, prep$dpars$mu)
})

# -----------------------------------------------------------------------
# ifelse() length-collapse regression -- see test-dlaplace1.R for the full
# explanation. dnorm2_lpmf_r's test (z >= mu) is built from z and mu, so a
# scalar z against vector mu/sigma is the direction that would collapse.
# -----------------------------------------------------------------------

test_that("dnorm2_lpmf_r does not collapse to length 1 for scalar z, vector mu/sigma", {
  z <- 2
  mu    <- c(-1, 0, 1, 2, 3)
  sigma <- c(1, 2, 3, 4, 5)
  out <- dnorm2_lpmf_r(z, mu, sigma)
  expect_length(out, length(mu))
  expect_equal(out, mapply(function(m, s) dnorm2_lpmf_r(z, m, s), mu, sigma))
})
