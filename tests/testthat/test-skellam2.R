# tests/testthat/test-skellam2.R
#
# Validates skellam2(), the asymmetric (free-mean) Skellam family, with
# the same rigor as skellam1's test-lpmf.R/test-lccdf.R: R-side formula
# vs skellam::dskellam/pskellam to near machine precision, Stan-side vs
# the same references once exposed, normalisation, branch-continuity at
# the production threshold seam, structural (non-rejection) enforcement
# of theta1/theta2 >= 0 confirmed via make_stancode(), and an end-to-end
# resp_trunc() recovery test with a genuinely nonzero mu_true.

# R-side log-PMF: same computation as log_lik_skellam2. sigma^2 = |mu| +
# sigmaexcess^2 -- see skellam2() in family.R for why this form (not the
# more obvious sigma = sqrt(mu^2 + sigmaexcess^2)) is required for
# theta1, theta2 >= 0.
r_lpmf2 <- function(k, mu, sigmaexcess) {
  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  z       <- 2 * sqrt(theta1 * theta2)
  log(besselI(z, abs(k), expon.scaled = TRUE)) + z - theta1 - theta2 + (k / 2) * log(theta1 / theta2)
}

# mu spans negative/zero/positive, including |mu| < 1 -- the regime
# where the rejected sigma = sqrt(mu^2 + sigmaexcess^2) construction
# fails (see family.R Details).
mu_vals          <- c(-20, -5, -0.5, 0, 0.5, 5, 20)
sigmaexcess_vals <- c(0.5, 1, 3, 10)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF matches dskellam across grid, including |mu| < 1", {
  k_vals <- -15L:15L

  for (mu in mu_vals) {
    for (se in sigmaexcess_vals) {
      sigmasq <- abs(mu) + se^2
      theta1  <- (sigmasq + mu) / 2
      theta2  <- (sigmasq - mu) / 2
      expect_true(theta1 >= 0 && theta2 >= 0,
                  label = paste0("mu=", mu, " se=", se, ": theta1=", theta1, " theta2=", theta2))

      r_vals   <- vapply(k_vals, r_lpmf2, numeric(1), mu = mu, sigmaexcess = se)
      ref_vals <- skellam::dskellam(k_vals, lambda1 = theta1, lambda2 = theta2, log = TRUE)
      expect_equal(r_vals, ref_vals, tolerance = 1e-6,
                   label = paste0("mu = ", mu, ", sigmaexcess = ", se))
    }
  }
})

test_that("R log-PMF normalises to 1 across grid", {
  for (mu in c(-5, -0.3, 0, 0.3, 5)) {
    for (se in c(0.5, 4)) {
      # Range scaled to the distribution's own spread (sigma), rather
      # than a fixed wide window -- summing to +-200 for a small-spread
      # case (e.g. se=0.5) evaluates besselI() at orders far beyond where
      # the PMF carries any mass, which triggers benign but noisy
      # "precision lost" warnings without changing the (negligible) sum.
      sigma  <- sqrt(abs(mu) + se^2)
      kmax   <- ceiling(abs(mu) + 20 * sigma + 20)
      total  <- suppressWarnings(sum(exp(r_lpmf2(-kmax:kmax, mu, se))))
      expect_equal(total, 1, tolerance = 1e-8,
                   label = paste0("mu = ", mu, ", sigmaexcess = ", se))
    }
  }
})

test_that("mu = 0 reduces exactly to skellam1's symmetric lpmf", {
  r_lpmf1 <- function(k, sigma) {
    mu <- sigma^2 / 2
    log(besselI(2 * mu, abs(k), expon.scaled = TRUE))
  }
  for (se in c(0.5, 2, 8)) {
    for (k in -5:5) {
      expect_equal(r_lpmf2(k, 0, se), r_lpmf1(k, se), tolerance = 1e-10,
                   label = paste0("se = ", se, ", k = ", k))
    }
  }
})

test_that("R log-PMF is numerically stable at large mu/sigmaexcess", {
  # Realistic-but-stressed range -- see the matching Stan lccdf stability
  # test below for the rationale (mu_hat in this project's real data
  # tops out around 30).
  for (mu in c(-50, -10, 0, 10, 50)) {
    for (sigma_target in c(10, 50, 100)) {
      se <- sqrt(max(0, sigma_target^2 - abs(mu)))
      for (k in round(mu) + c(-5L, 0L, 5L)) {
        val <- r_lpmf2(k, mu, se)
        expect_false(is.nan(val),      label = paste0("NaN at mu=", mu, " sigma=", sigma_target, " k=", k))
        expect_false(is.infinite(val), label = paste0("Inf at mu=", mu, " sigma=", sigma_target, " k=", k))
      }
    }
  }
})

log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

logS_bruteforce2 <- function(y, mu, sigmaexcess, K = NULL) {
  sigmasq <- abs(mu) + sigmaexcess^2
  if (is.null(K)) K <- ceiling(abs(mu) + sigmasq + 40 * sqrt(sigmasq) + 50)
  # suppressWarnings: besselI() flags "precision lost" whenever the order
  # is far beyond where the argument gives it any real mass (e.g. a small
  # se combined with a deliberately large y in this grid) -- the returned
  # value is still a correctly-negligible contribution to the sum, this
  # is just R being conservative about a term that doesn't matter.
  suppressWarnings(log_sum_exp(r_lpmf2((-y):K, mu, sigmaexcess)))
}

logS_pkg2 <- function(y, mu, sigmaexcess) {
  sigmasq <- abs(mu) + sigmaexcess^2
  theta1  <- (sigmasq + mu) / 2
  theta2  <- (sigmasq - mu) / 2
  log1p(-skellam::pskellam(-y - 1, lambda1 = theta1, lambda2 = theta2))
}

trunc_grid <- expand.grid(
  mu          = c(-20, -5, 0, 5, 20),
  sigmaexcess = c(0.5, 1, 5, 10),
  y           = c(0, 1, 3, 10, 20)
)

test_that("R brute-force log-CCDF agrees with skellam::pskellam across grid", {
  diffs <- mapply(
    function(mu, se, y) logS_pkg2(y, mu, se) - logS_bruteforce2(y, mu, se),
    trunc_grid$mu, trunc_grid$sigmaexcess, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-7)
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lpmf_stan_code <- paste0("functions {\n", skellam2_stan_funs, "}\nmodel {}\n")

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

test_that("Stan log-PMF matches R log-PMF across grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  k_vals <- -10L:10L

  for (mu in mu_vals) {
    for (se in sigmaexcess_vals) {
      stan_vals <- vapply(k_vals, function(k) skellam2_lpmf(k, mu, se), numeric(1))
      r_vals    <- vapply(k_vals, r_lpmf2, numeric(1), mu = mu, sigmaexcess = se)
      expect_equal(stan_vals, r_vals, tolerance = 1e-6,
                   label = paste0("mu = ", mu, ", sigmaexcess = ", se))
    }
  }
})

lccdf_stan_block2 <- function(threshold = 100) {
  paste0("functions {\n", skellam2_lccdf_stan(threshold), "}\nmodel {}\n")
}

lccdf_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm_lccdf <- rstan::stan_model(model_code = lccdf_stan_block2(100))
      rstan::expose_stan_functions(sm_lccdf)
    })
    lccdf_ready <- TRUE
  }, error = function(e) NULL)
}

# Grid chosen so (theta1 + theta2) / 2 stays comfortably under the
# default threshold (100) for every combination -- max here is
# (20 + 100) / 2 = 60 -- so every point exercises the exact branch
# unambiguously (see skellam1's test-lccdf.R for why landing exactly on
# the threshold is avoided: floating-point round-trips near a strict
# inequality boundary are inherently brittle, not a logic bug).
test_that("Stan skellam2_lccdf (exact branch) matches R brute-force reference", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  diffs <- mapply(
    function(mu, se, y) skellam2_lccdf(as.integer(-y - 1), mu, se) - logS_bruteforce2(y, mu, se),
    trunc_grid$mu, trunc_grid$sigmaexcess, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam2_lccdf matches package pskellam reference", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  diffs <- mapply(
    function(mu, se, y) skellam2_lccdf(as.integer(-y - 1), mu, se) - logS_pkg2(y, mu, se),
    trunc_grid$mu, trunc_grid$sigmaexcess, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam2_lccdf is numerically stable at large sigma (normal-approx branch)", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  # Realistic-but-stressed range: in the count data these families are built
  # for, per-group mu_hat tops out around 30, so mu in [-50, 50] and
  # sigma in [0, 100] is comfortably beyond anything this model will see
  # without being an arbitrary, disconnected extreme. y is chosen near
  # each mu (where real data would actually land), not deliberately far
  # from it.
  for (mu in c(-50, -10, 0, 10, 50)) {
    for (sigma_target in c(10, 50, 100)) {
      se <- sqrt(max(0, sigma_target^2 - abs(mu)))
      for (y in round(mu) + c(-5L, 0L, 5L)) {
        val <- skellam2_lccdf(y, mu, se)
        expect_false(is.nan(val),      label = paste0("NaN at mu=", mu, " sigma=", sigma_target, " y=", y))
        expect_false(is.infinite(val), label = paste0("Inf at mu=", mu, " sigma=", sigma_target, " y=", y))
      }
    }
  }
})

test_that("exact and normal-approx branches agree closely AT the production cutover seam", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  exact_always_env  <- new.env()
  approx_always_env <- new.env()

  exact_ready <- tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block2(1e15))
      rstan::expose_stan_functions(sm, env = exact_always_env)
    })
    TRUE
  }, error = function(e) FALSE)

  approx_ready <- tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block2(-1))
      rstan::expose_stan_functions(sm, env = approx_always_env)
    })
    TRUE
  }, error = function(e) FALSE)

  skip_if_not(exact_ready && approx_ready, "Stan compilation of seam-test variants failed")

  threshold      <- 100
  mu_skellam_seam <- threshold + c(-10, -1, 0, 1, 10)
  mu_fixed_vals  <- c(-10, 0, 10)
  y_vals         <- c(0L, 3L, 10L)

  seam_grid <- expand.grid(mu = mu_fixed_vals, mu_skellam = mu_skellam_seam, y = y_vals)
  # Solve sigmaexcess so that (theta1+theta2)/2 = (|mu|+sigmaexcess^2)/2
  # lands exactly on each seam point, for each fixed mu.
  seam_grid$sigmaexcess <- sqrt(2 * seam_grid$mu_skellam - abs(seam_grid$mu))

  diffs <- mapply(
    function(mu, se, y) {
      exact_val  <- exact_always_env$skellam2_lccdf(y, mu, se)
      approx_val <- approx_always_env$skellam2_lccdf(y, mu, se)
      abs(exact_val - approx_val)
    },
    seam_grid$mu, seam_grid$sigmaexcess, seam_grid$y
  )

  expect_lt(max(diffs), 0.01,
            label = paste0("max diff near seam = ", round(max(diffs), 5)))
})

# -----------------------------------------------------------------------
# Structural-enforcement check: theta1, theta2 >= 0 by construction
# -----------------------------------------------------------------------

test_that("make_stancode() shows theta1/theta2 enforced structurally, not via rejection", {
  skip_if_not_installed("brms")

  fake_dat <- data.frame(delta = sample(-5:5, 20, replace = TRUE), x = rnorm(20))
  code <- brms::make_stancode(
    brms::bf(delta ~ x, sigmaexcess ~ x),
    family   = skellam2(),
    stanvars = skellam2_stanvars(),
    data     = fake_dat
  )

  # No reject() statement anywhere -- theta1/theta2 positivity comes from
  # the sigma^2 = |mu| + sigmaexcess^2 construction inside skellam2_lpmf,
  # not from sampler-level rejection of invalid draws.
  expect_false(grepl("reject\\(", code))
  # The deterministic construction itself is present in the generated
  # functions block.
  expect_true(grepl("abs(mu) + square(sigmaexcess)", code, fixed = TRUE))
})

# -----------------------------------------------------------------------
# End-to-end smoke test: resp_trunc() recovery with a nonzero mu_true
# -----------------------------------------------------------------------

test_that("skellam2_stanvars() + skellam2_lccdf_stanvars() recovers nonzero mu_true and sigma_true under resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  set.seed(2)
  n          <- 60
  mu_true    <- 3
  sigma_true <- 5
  theta1_true <- (sigma_true^2 + mu_true) / 2
  theta2_true <- (sigma_true^2 - mu_true) / 2
  y2       <- sample(0:10, n, replace = TRUE)

  draw_trunc <- function(y) {
    repeat {
      d <- skellam::rskellam(1, lambda1 = theta1_true, lambda2 = theta2_true)
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y2, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y2)

  # Same documented default-prior fragility as skellam1 (see
  # test-lccdf.R) -- applying weakly-informative priors on both
  # intercepts up front rather than discovering the same failure mode again.
  sane_prior <- brms::prior(normal(0, 5), class = "Intercept") +
    brms::prior(normal(1, 1.5), class = "Intercept", dpar = "sigmaexcess")

  suppressMessages({
    fit <- brms::brm(
      brms::bf(delta | trunc(lb = neg_bound) ~ 1, sigmaexcess ~ 1),
      family   = skellam2(),
      stanvars = skellam2_stanvars() + skellam2_lccdf_stanvars(),
      data     = dat,
      prior    = sane_prior,
      backend  = "cmdstanr",
      chains   = 2,
      iter     = 1000,
      warmup   = 500,
      seed     = 2,
      refresh  = 0
    )
  })

  expect_s3_class(fit, "brmsfit")

  draws <- as.data.frame(fit)

  mu_hat    <- draws[["b_Intercept"]]
  mu_q      <- quantile(mu_hat, c(0.025, 0.975))
  expect_true(
    mu_true >= mu_q[[1]] && mu_true <= mu_q[[2]],
    label = paste0("mu_true = ", mu_true, ", 95% CI: [", round(mu_q[[1]], 3), ", ", round(mu_q[[2]], 3), "]")
  )

  sigmaexcess_hat <- exp(draws[["b_sigmaexcess_Intercept"]])
  sigma_hat       <- sqrt(abs(mu_hat) + sigmaexcess_hat^2)
  sigma_q         <- quantile(sigma_hat, c(0.025, 0.975))
  expect_true(
    sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    label = paste0("sigma_true = ", sigma_true, ", 95% CI: [", round(sigma_q[[1]], 3), ", ", round(sigma_q[[2]], 3), "]")
  )

  # Extend the existing fit (no new brm() call) to also check that
  # posterior_predict()/posterior_epred() honour resp_trunc()'s bounds --
  # confirms the R-side truncation fix works through the full brms
  # dispatch path, not just the synthetic-prep unit tests in this file.
  pp <- brms::posterior_predict(fit)
  expect_true(all(sweep(pp, 2, dat$neg_bound, `>=`)),
              label = "posterior_predict draws below their row's trunc(lb=) bound")

  # brms::posterior_epred() itself cannot be used here -- see the matching
  # comment in test-lccdf.R for the confirmed brms-level reason
  # (posterior_epred.brmsprep() unconditionally routes truncated fits to
  # posterior_epred_trunc(), which has no generic handler for ANY custom
  # family and always errors, independent of this package's own fix).
  prep <- brms::prepare_predictions(fit)
  ep   <- posterior_epred_skellam2(prep)
  expect_true(all(is.finite(ep)))
  expect_true(all(sweep(ep, 2, dat$neg_bound, `>=`)),
              label = "posterior_epred below the row's trunc(lb=) bound")
})

# -----------------------------------------------------------------------
# posterior_predict / posterior_epred truncation correctness (synthetic prep)
# -----------------------------------------------------------------------

test_that("posterior_predict_skellam2 respects lb (repro of confirmed bug)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 1500, ncol = 1),
                 sigmaexcess = matrix(2, nrow = 1500, ncol = 1)),
    Y = 0, lb = 0
  )
  set.seed(123)
  draws <- posterior_predict_skellam2(1, prep)
  expect_true(all(draws >= 0), label = paste0("min draw = ", min(draws)))
})

test_that("posterior_predict_skellam2 without lb/ub matches untruncated behaviour (fast path)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 3000, ncol = 1),
                 sigmaexcess = matrix(2, nrow = 3000, ncol = 1)),
    Y = 0
  )
  set.seed(1)
  draws <- posterior_predict_skellam2(1, prep)
  expect_true(is.numeric(draws) && length(draws) == 3000)
  expect_true(min(draws) < 0)
})

test_that("posterior_predict_skellam2 draws match the truncated PMF distributionally", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 4000, ncol = 1),
                 sigmaexcess = matrix(2, nrow = 4000, ncol = 1)),
    Y = 0, lb = 0, ub = 20
  )
  set.seed(42)
  draws <- posterior_predict_skellam2(1, prep)
  expect_true(all(draws >= 0 & draws <= 20))

  support <- 0:20
  probs <- exp(skellam2_lpmf_r(support, 4, 2))
  probs <- probs / sum(probs)
  emp <- as.numeric(table(factor(draws, levels = support))) / length(draws)
  expect_lt(max(abs(emp - probs)), 0.03)
})

test_that("posterior_epred_skellam2 differs from untruncated mu when lb is tight", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 5, ncol = 1),
                 sigmaexcess = matrix(2, nrow = 5, ncol = 1)),
    Y = 0, lb = 0
  )
  epred <- posterior_epred_skellam2(prep)
  expect_equal(dim(epred), c(5, 1))
  expect_true(all(epred[, 1] > 4), label = paste0("epred = ", epred[1, 1]))
})

test_that("posterior_epred_skellam2 matches brute-force truncated mean", {
  mu_val <- 4; sigmaexcess_val <- 2; lb_val <- 0
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(mu_val, nrow = 1, ncol = 1),
                 sigmaexcess = matrix(sigmaexcess_val, nrow = 1, ncol = 1)),
    Y = 0, lb = lb_val
  )
  epred <- posterior_epred_skellam2(prep)

  support <- lb_val:80
  probs <- exp(skellam2_lpmf_r(support, mu_val, sigmaexcess_val))
  probs <- probs / sum(probs)
  brute_force_mean <- sum(support * probs)

  expect_equal(epred[1, 1], brute_force_mean, tolerance = 1e-6)
})

test_that("posterior_epred_skellam2 leaves untruncated observations exactly at mu (no regression)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(c(4, -3), nrow = 2, ncol = 2),
                 sigmaexcess = matrix(2, nrow = 2, ncol = 2)),
    Y = c(0, 0)
  )
  epred <- posterior_epred_skellam2(prep)
  expect_equal(epred, prep$dpars$mu)
})
