# tests/testthat/test-lccdf.R
#
# R-side and Stan-side validation of skellam1_lccdf: agreement against
# skellam::pskellam() and a
# brute-force log_sum_exp tail-sum, both on the R side and, once exposed
# via rstan::expose_stan_functions(), on the Stan side too.
#
# Reference helpers (lccdf_r_lpmf, logS_bruteforce, logS_pkg) are kept on
# the mu_skellam scale -- the underlying Skellam rate -- since that's the
# natural scale for the Bessel-sum/CDF maths being checked. The Stan-side
# skellam1_lccdf() now takes sigma (sigma = sqrt(2 * mu_skellam)) as its
# second argument, so Stan-side tests convert at the call site.

# R-side log-PMF (same computation as log_lik_skellam1 / r_lpmf in
# test-lpmf.R; duplicated locally so this file has no cross-file
# dependency on test execution order).
lccdf_r_lpmf <- function(k, mu) log(besselI(2 * mu, abs(k), expon.scaled = TRUE))

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# log S(y, mu) = log P(delta >= -y), the truncated normalising constant.
# lccdf(k, mu) = log P(delta > k), so logS(y, mu) = lccdf(-y - 1, mu).
logS_bruteforce <- function(y, mu, K = NULL) {
  if (is.null(K)) K <- ceiling(mu + 40 * sqrt(mu) + 50)
  log_sum_exp(lccdf_r_lpmf((-y):K, mu))
}

logS_pkg <- function(y, mu) {
  log1p(-skellam::pskellam(-y - 1, lambda1 = mu, lambda2 = mu))
}

mu_grid    <- c(0.2, 0.5, 1, 2, 5, 10, 25, 50, 100)
y_grid     <- c(0, 1, 2, 3, 5, 10, 20, 30)
grid       <- expand.grid(mu = mu_grid, y = y_grid)
sigma_grid <- sqrt(2 * mu_grid)

# Separate grid for Stan-side exact-branch tests: 100 is replaced by 99.
# skellam1_lccdf takes sigma and derives mu_skellam = sigma^2 / 2
# internally; at mu_skellam exactly 100 (the default normal_approx_threshold),
# the sqrt(2*mu) -> square()/2 round-trip lands at 100.0000000000000142 (a
# ~1.4e-14 floating-point perturbation), which flips the strict `>` branch
# comparison to the normal-approx side even though the logical input is
# the boundary value itself. This is a genuine, confirmed floating-point
# edge effect of deriving mu_skellam from sigma rather than supplying it
# directly (verified: exact and approx already agree to ~1e-4 at this
# mu_skellam regardless -- see the seam-continuity test below -- so this
# is immaterial in practice, but it does mean a precision test landing
# exactly on the threshold value is testing branch-selection brittleness,
# not the lpmf/lccdf maths). Testing at 99 keeps this grid's intent (probe
# close to the threshold while staying unambiguously on the exact-branch
# side) without that brittleness.
mu_grid_stan <- c(0.2, 0.5, 1, 2, 5, 10, 25, 50, 99)
grid_stan    <- expand.grid(mu = mu_grid_stan, y = y_grid)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force log-CCDF agrees with skellam::pskellam across grid", {
  diffs <- mapply(
    function(mu, y) logS_pkg(y, mu) - logS_bruteforce(y, mu),
    grid$mu, grid$y
  )
  expect_lt(max(abs(diffs)), 1e-8)
})

test_that("leakage at y = 0 matches the documented structural finding", {
  # At y = 0, leakage = 1 - S(0, mu) should be at least ~15% for any mu --
  # the key structural finding motivating truncation in the first place.
  leakage <- vapply(mu_grid, function(mu) 1 - exp(logS_bruteforce(0, mu)), numeric(1))
  expect_true(all(leakage > 0.14),
              label = paste0("min leakage = ", round(min(leakage), 3)))
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lccdf_stan_block <- function(threshold = 100) {
  paste0("functions {\n", skellam1_lccdf_stan(threshold), "}\nmodel {}\n")
}

stan_ready <- FALSE
if (requireNamespace("rstan", quietly = TRUE)) {
  tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block(100))
      rstan::expose_stan_functions(sm)
    })
    stan_ready <- TRUE
  }, error = function(e) NULL)
}

test_that("Stan skellam1_lccdf (exact branch) matches R brute-force reference", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  # Grid is entirely mu < 100, so every point uses the exact (Bessel-sum)
  # branch under the default normal_approx_threshold (on the mu_skellam
  # scale derived from sigma inside skellam1_lccdf) unambiguously -- see
  # mu_grid_stan above for why 100 itself is excluded here.
  diffs <- mapply(
    function(mu, y) skellam1_lccdf(as.integer(-y - 1), sqrt(2 * mu)) - logS_bruteforce(y, mu),
    grid_stan$mu, grid_stan$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam1_lccdf matches package pskellam reference", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  diffs <- mapply(
    function(mu, y) skellam1_lccdf(as.integer(-y - 1), sqrt(2 * mu)) - logS_pkg(y, mu),
    grid_stan$mu, grid_stan$y
  )
  expect_lt(max(abs(diffs)), 1e-6)
})

test_that("Stan skellam1_lccdf is numerically stable at large sigma (normal-approx branch)", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  for (mu in c(500, 5000)) {
    sigma <- sqrt(2 * mu)
    for (y in c(-5L, 0L, 5L)) {
      val <- skellam1_lccdf(y, sigma)
      expect_false(is.nan(val),      label = paste0("NaN at sigma=", sigma, " y=", y))
      expect_false(is.infinite(val), label = paste0("Inf at sigma=", sigma, " y=", y))
    }
  }
})

test_that("normal_approx_threshold parameter actually moves the cutover point", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  low_thresh_env <- new.env()
  ok <- tryCatch({
    suppressMessages({
      sm_low <- rstan::stan_model(model_code = lccdf_stan_block(5))
      rstan::expose_stan_functions(sm_low, env = low_thresh_env)
    })
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "Stan compilation of custom-threshold variant failed")

  mu    <- 10  # > 5 (custom threshold) but <= 100 (default threshold), on mu_skellam scale
  sigma <- sqrt(2 * mu)
  y     <- 3L

  exact_val  <- skellam1_lccdf(y, sigma)                  # default threshold: exact branch
  approx_val <- low_thresh_env$skellam1_lccdf(y, sigma)   # threshold = 5: normal-approx branch

  expect_false(is.nan(approx_val))
  expect_false(is.infinite(approx_val))
  # Different branches at the same (y, sigma) should give visibly different
  # results -- confirms the parameter actually moves the cutover, not
  # just that the generated code still parses.
  expect_gt(abs(exact_val - approx_val), 1e-4)
})

test_that("exact and normal-approx branches agree closely AT the production cutover seam", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  # The "different threshold values produce different results" test above
  # confirms the cutover point moves, but not that the two branches agree
  # with each other near the seam they share. Force each branch
  # unconditionally (threshold = 1e15 always takes the exact branch for
  # any mu_skellam tested here; threshold = -1 always takes the
  # normal-approx branch, since mu_skellam >= 0 > -1 always) and compare
  # them directly at mu_skellam = 100 +/- {0, 1, 10} -- the seam around
  # the package's actual default threshold of 100.
  exact_always_env  <- new.env()
  approx_always_env <- new.env()

  exact_ready <- tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block(1e15))
      rstan::expose_stan_functions(sm, env = exact_always_env)
    })
    TRUE
  }, error = function(e) FALSE)

  approx_ready <- tryCatch({
    suppressMessages({
      sm <- rstan::stan_model(model_code = lccdf_stan_block(-1))
      rstan::expose_stan_functions(sm, env = approx_always_env)
    })
    TRUE
  }, error = function(e) FALSE)

  skip_if_not(exact_ready && approx_ready, "Stan compilation of seam-test variants failed")

  threshold  <- 100
  mu_seam    <- threshold + c(-10, -1, 0, 1, 10)
  sigma_seam <- sqrt(2 * mu_seam)
  y_vals     <- c(0L, 3L, 10L)

  seam_grid <- expand.grid(sigma = sigma_seam, y = y_vals)
  diffs <- mapply(
    function(sigma, y) {
      exact_val  <- exact_always_env$skellam1_lccdf(y, sigma)
      approx_val <- approx_always_env$skellam1_lccdf(y, sigma)
      abs(exact_val - approx_val)
    },
    seam_grid$sigma, seam_grid$y
  )

  # Empirically, max abs diff near this seam is ~8e-4 (y=10, mu=90); 0.01
  # is a generous but still meaningful bound -- a real formula bug (sign
  # error, wrong variance, off-by-one in the continuity correction) would
  # blow well past it.
  expect_lt(max(diffs), 0.01,
            label = paste0("max diff near seam = ", round(max(diffs), 5)))
})

# -----------------------------------------------------------------------
# End-to-end smoke test: the resp_trunc() usage example from the docs
# -----------------------------------------------------------------------

test_that("skellam1_stanvars() + skellam1_lccdf_stanvars() works with resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  set.seed(1)
  n          <- 40
  mu_true    <- 5
  sigma_true <- sqrt(2 * mu_true)
  y2       <- sample(0:8, n, replace = TRUE)

  # Rejection-sample from the truncated Skellam(mu_true, mu_true) subject
  # to delta >= -y2, to build a dataset whose generative process matches
  # what resp_trunc() should correct for.
  draw_trunc <- function(y) {
    repeat {
      d <- skellam::rskellam(1, lambda1 = mu_true, lambda2 = mu_true)
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y2, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y2)

  # brms's default Intercept prior (student_t(3, 0, 2.5)) is known to be
  # wide enough that this custom Bessel-based likelihood can occasionally
  # wander to a nonsensical log(sigma) region for an unlucky seed -- hit
  # here for seed=1 under the sigma-reparameterisation (165 divergences,
  # Rhat 1.84 without it). Applying that same documented fix.
  sane_prior <- brms::prior(normal(1, 1.5), class = "Intercept")

  suppressMessages({
    fit <- brms::brm(
      brms::bf(delta | trunc(lb = neg_bound) ~ 1),
      family   = skellam1(),
      stanvars = skellam1_stanvars() + skellam1_lccdf_stanvars(),
      data     = dat,
      prior    = sane_prior,
      backend  = "cmdstanr",
      chains   = 2,
      iter     = 800,
      warmup   = 400,
      seed     = 1,
      refresh  = 0
    )
  })

  expect_s3_class(fit, "brmsfit")

  # CI-coverage standard (matching the package's own lpmf recovery tests
  # in test-recovery.R), not just a finiteness check: sigma_true should
  # fall within the 95% posterior CrI for sigma = exp(Intercept) (log link).
  draws     <- as.data.frame(fit)
  sigma_hat <- exp(draws[["b_Intercept"]])
  sigma_q   <- quantile(sigma_hat, c(0.025, 0.975))
  expect_true(
    sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    label = paste0("sigma_true = ", round(sigma_true, 3),
                   ", 95% CI: [", round(sigma_q[[1]], 3),
                   ", ", round(sigma_q[[2]], 3), "]")
  )

  # Extend the existing fit (no new brm() call) to also check that
  # posterior_predict()/posterior_epred() honour resp_trunc()'s bounds --
  # confirms the R-side truncation fix (see R/truncation.R and
  # posterior_predict_skellam1/posterior_epred_skellam1 in R/family.R)
  # works through the full brms dispatch path, not just the synthetic-prep
  # unit tests in test-skellam1.R.
  pp <- brms::posterior_predict(fit)
  expect_true(all(sweep(pp, 2, dat$neg_bound, `>=`)),
              label = "posterior_predict draws below their row's trunc(lb=) bound")

  # brms::posterior_epred() itself cannot be used here: brms's generic
  # posterior_epred.brmsprep() checks is_trunc() *before* considering
  # family type and unconditionally routes truncated fits to
  # posterior_epred_trunc(), which looks for a function named
  # "posterior_epred_trunc_custom" inside brms's OWN namespace (not this
  # package's) -- confirmed via brms:::posterior_epred_trunc's source.
  # That function doesn't exist for any custom family, so
  # brms::posterior_epred() always errors ("not yet implemented for
  # truncated 'custom' models") on a truncated custom-family fit,
  # regardless of whether the family's own posterior_epred_<name>() is
  # correct. Calling posterior_epred_skellam1() directly on a real
  # prepare_predictions() object sidesteps that brms-level gate while
  # still exercising the real fit end-to-end.
  prep <- brms::prepare_predictions(fit)
  ep   <- posterior_epred_skellam1(prep)
  expect_true(all(is.finite(ep)))
  expect_true(all(sweep(ep, 2, dat$neg_bound, `>=`)),
              label = "posterior_epred below the row's trunc(lb=) bound")
})
