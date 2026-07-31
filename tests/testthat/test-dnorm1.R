# tests/testthat/test-dnorm1.R
#
# Validates dnorm1(), the discrete normal family (location fixed at 0,
# free scale), with the same rigor as dlaplace1/dlaplace2's test suites.
#
# Reference note: unlike dlaplace1 (where extraDistr::ddlaplace() turned
# out to be a *different* distribution entirely, forcing a hand-derived
# reference), base R's pnorm() is directly, correctly parameterised for
# this family -- pnorm(x, sd = sigma) is exactly the continuous Normal(0,
# sigma) CDF this family discretises. No reference mismatch risk here.
#
# Numerical-stability note, confirmed and load-bearing for this entire
# file: a naive log(pnorm(z+0.5,sd=sigma) - pnorm(z-0.5,sd=sigma)) (or
# its log-space cousin using pnorm(..., log.p = TRUE) on both terms)
# catastrophically cancels once z is far enough into the positive tail
# that both pnorm() calls round to the same double -- confirmed to occur
# at only ~10 SDs out for sigma=1, *inside* this package's own
# realistic-but-stressed test range for the other families (sigma up to
# 100, k within 10 SDs of the centre). The reference below fixes this by
# differencing two *survival* values (pnorm(..., lower.tail=FALSE)) when
# z is on the far side of the mean, mirroring dnorm1_lpmf's own fix in
# stanfunctions.R (see there for the full derivation) -- this is not
# optional defensive coding, it is required for this reference to be
# trustworthy across the grids used below.
#
# Separately, and importantly: Stan's *built-in* normal_lccdf is not a
# safe building block for the survival side of that fix, despite being
# the obvious choice. This is a documented Stan limitation, not a guess:
# the Stan Functions Reference states normal_lccdf underflows to -inf
# for (y-mu)/sigma above ~8.25, and stan-dev/math issue #1985
# (https://github.com/stan-dev/math/issues/1985) confirms normal_lccdf
# (unlike normal_lcdf, which received it) was never updated with the
# more accurate Mills-ratio approximation, so it still hits this floor.
# Confirmed directly during this family's development: normal_lccdf(9.5
# | 0, 1) returns -inf here vs the true -48.3. dnorm1_lpmf/dnorm1_lccdf
# in stanfunctions.R use the exact closed form P(Z>x) =
# 0.5*erfc(x/(sigma*sqrt(2))) instead -- confirmed to match R's pnorm()
# to high precision out to 30+ SDs, since erfc() does not share
# normal_lccdf's accuracy gap. The same documented issue, independently
# rediscovered, motivated reinstating an equivalent erfc-based fix in
# skellam1_lccdf/skellam2_lccdf's normal-approximation branch (see
# stanfunctions.R) after an earlier, unrelated investigation had reverted
# it as seemingly unnecessary.

log_normal_cdf   <- function(x, sigma) pnorm(x, sd = sigma, log.p = TRUE)
log_normal_lccdf <- function(x, sigma) pnorm(x, sd = sigma, lower.tail = FALSE, log.p = TRUE)
log_diff_exp_r   <- function(a, bb) a + log1p(-exp(bb - a))

# R-side log-PMF: stable log-space reference, branched on the sign of z
# for the reason documented above.
r_lpmf_dnorm1 <- function(z, sigma) {
  ifelse(
    z >= 0,
    log_diff_exp_r(log_normal_lccdf(z - 0.5, sigma), log_normal_lccdf(z + 0.5, sigma)),
    log_diff_exp_r(log_normal_cdf(z + 0.5, sigma), log_normal_cdf(z - 0.5, sigma))
  )
}

sigma_vals <- c(0.2, 0.5, 1, 2, 5, 10, 25, 50, 100)

# -----------------------------------------------------------------------
# R-side tests — no Stan compilation required
# -----------------------------------------------------------------------

test_that("R log-PMF normalises to 1 across grid", {
  for (sigma in sigma_vals) {
    # kmax anchored to ~10 SDs (not dlaplace1's 40*sigma+50): confirmed
    # empirically that the discrete normal's much thinner (Gaussian, not
    # exponential) tail needs far fewer SDs to capture essentially all
    # the mass -- at sigma=100, 8 SDs (kmax=800) already gives a total
    # mass of 1 to 12 decimal places; 10 SDs is a comfortable margin
    # beyond that, not a value blindly carried over from dlaplace1.
    kmax  <- ceiling(10 * sigma + 20)
    total <- sum(exp(r_lpmf_dnorm1(-kmax:kmax, sigma)))
    expect_equal(total, 1, tolerance = 1e-8, label = paste0("sigma = ", sigma))
  }
})

test_that("R log-PMF is numerically stable, realistic-but-stressed range", {
  # Mirrors the realistic-but-stressed range agreed for skellam2/
  # dlaplace1/dlaplace2 (sigma up to 100, k within 10 SDs) -- this range
  # is anchored to this project's actual data scale (real per-taxon
  # mu_hat / sigma_hat tops out well under 50), not to dlaplace1's own
  # convention, but happens to use the same numbers since both are
  # expressing "comfortably beyond real data, not an arbitrary extreme."
  for (sigma in c(0.2, 1, 10, 100)) {
    for (k in round(sigma * c(-10, 0, 10))) {
      val <- r_lpmf_dnorm1(k, sigma)
      expect_false(is.nan(val),      label = paste0("NaN at sigma=", sigma, " k=", k))
      expect_false(is.infinite(val), label = paste0("Inf at sigma=", sigma, " k=", k))
    }
  }
})

test_that("R log-PMF stays finite well beyond the realistic range too (branch fix confirmed)", {
  # Without the z>=0 branch in r_lpmf_dnorm1, this grid reproduces the
  # exact failure this family's Stan implementation also had to be
  # fixed for: -inf by 10 SDs out. Confirms the fix holds much further
  # (the normal's tail is thin enough that even 30+ SDs stays finite
  # here, unlike dlaplace1's reference, which has its own, much earlier,
  # genuine double-precision floor).
  for (sigma in c(0.2, 1, 10, 100)) {
    for (k_mult in c(15, 20, 25, 30)) {
      val <- r_lpmf_dnorm1(round(sigma * k_mult), sigma)
      expect_false(is.infinite(val), label = paste0("sigma=", sigma, " k_mult=", k_mult))
    }
  }
})

log_sum_exp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# Brute-force log P(Z > y): sums from y+1, matching dnorm1_lccdf's direct
# log P(Z>y) definition (the same "not the skellam-style S(y)=P(>=-y)
# quantity" caveat noted in test-dlaplace1.R applies here too).
logS_bruteforce_dnorm1 <- function(y, sigma, K = NULL) {
  if (is.null(K)) K <- y + ceiling(10 * sigma + 20)
  log_sum_exp(r_lpmf_dnorm1((y + 1):K, sigma))
}

# Closed-form log P(Z > y) via the same stable survival function the
# reference above uses for its z>=0 branch.
logS_closed_dnorm1 <- function(y, sigma) log_normal_lccdf(y + 0.5, sigma)

# y as a multiple of sigma (SD units) -- pushed out to 10 SDs, well
# beyond what the brute-force / closed-form comparison needs (8 SDs
# already captures the mass to 12 decimal places; 10 leaves margin).
y_multiples <- c(0, 1, 2, 3, 5, 10)
trunc_grid  <- expand.grid(sigma = sigma_vals, y_mult = y_multiples)
trunc_grid$y <- round(trunc_grid$sigma * trunc_grid$y_mult)

test_that("R brute-force log-CCDF agrees with the closed form across grid", {
  diffs <- mapply(
    function(sigma, y) logS_closed_dnorm1(y, sigma) - logS_bruteforce_dnorm1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs)), 1e-8)
})

# -----------------------------------------------------------------------
# Stan tests — require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

lpmf_stan_code <- paste0("functions {\n", dnorm1_stan_funs, "}\nmodel {}\n")

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
  for (sigma in sigma_vals) {
    k_vals    <- unique(round(sigma * (-10:10)))
    stan_vals <- vapply(k_vals, function(k) dnorm1_lpmf(as.integer(k), sigma), numeric(1))
    r_vals    <- vapply(k_vals, r_lpmf_dnorm1, numeric(1), sigma = sigma)
    expect_equal(stan_vals, r_vals, tolerance = 1e-6, label = paste0("sigma = ", sigma))
  }
})

test_that("Stan log-PMF stays finite well beyond 10 SDs (erfc fix confirmed, not just the R reference)", {
  # This is the test that would have failed before the erfc fix: the
  # branch alone (differencing survival values instead of CDF values)
  # is not sufficient if the survival values themselves come from the
  # documented-broken normal_lccdf -- see the file header citation.
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")
  for (sigma in c(1, 10, 100)) {
    for (k_mult in c(9, 15, 20, 25)) {
      val <- dnorm1_lpmf(as.integer(round(sigma * k_mult)), sigma)
      expect_false(is.infinite(val), label = paste0("sigma=", sigma, " k_mult=", k_mult))
    }
  }
})

lccdf_stan_code <- paste0("functions {\n", dnorm1_lccdf_stan, "}\nmodel {}\n")

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

test_that("Stan dnorm1_lccdf matches R closed-form / brute-force references, including past 8.25 SDs", {
  skip_if_not(lccdf_ready, "rstan unavailable or Stan compilation failed")
  diffs_closed <- mapply(
    function(sigma, y) dnorm1_lccdf(as.integer(y), sigma) - logS_closed_dnorm1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_closed)), 1e-6)

  diffs_brute <- mapply(
    function(sigma, y) dnorm1_lccdf(as.integer(y), sigma) - logS_bruteforce_dnorm1(y, sigma),
    trunc_grid$sigma, trunc_grid$y
  )
  expect_lt(max(abs(diffs_brute)), 1e-6)
})

test_that("Stan dnorm1_lpmf and dnorm1_lccdf are numerically stable, realistic-but-stressed range", {
  skip_if_not(stan_ready && lccdf_ready, "rstan unavailable or Stan compilation failed")
  for (sigma in c(0.2, 1, 10, 100)) {
    for (k in round(sigma * c(-10, 0, 10))) {
      lpmf_val  <- dnorm1_lpmf(as.integer(k), sigma)
      lccdf_val <- dnorm1_lccdf(as.integer(k), sigma)
      expect_false(is.nan(lpmf_val) || is.infinite(lpmf_val),
                   label = paste0("lpmf at sigma=", sigma, " k=", k))
      expect_false(is.nan(lccdf_val) || is.infinite(lccdf_val),
                   label = paste0("lccdf at sigma=", sigma, " k=", k))
    }
  }
})

# -----------------------------------------------------------------------
# End-to-end smoke test: resp_trunc() recovery
# -----------------------------------------------------------------------

test_that("dnorm1_stanvars() + dnorm1_lccdf_stanvars() recovers sigma_true under resp_trunc()", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  set.seed(7)
  n          <- 50
  sigma_true <- 6
  y2       <- sample(0:8, n, replace = TRUE)

  draw_trunc <- function(y) {
    repeat {
      d <- round(rnorm(1, mean = 0, sd = sigma_true))
      if (d >= -y) return(d)
    }
  }
  delta <- vapply(y2, draw_trunc, numeric(1))
  dat   <- data.frame(delta = delta, neg_bound = -y2)

  sane_prior <- brms::prior(normal(1, 1.5), class = "Intercept")

  suppressMessages({
    fit <- brms::brm(
      brms::bf(delta | trunc(lb = neg_bound) ~ 1),
      family   = dnorm1(),
      stanvars = dnorm1_stanvars() + dnorm1_lccdf_stanvars(),
      data     = dat,
      prior    = sane_prior,
      backend  = "cmdstanr",
      chains   = 2,
      iter     = 800,
      warmup   = 400,
      seed     = 7,
      refresh  = 0
    )
  })

  expect_s3_class(fit, "brmsfit")

  draws     <- as.data.frame(fit)
  sigma_hat <- exp(draws[["b_Intercept"]])
  sigma_q   <- quantile(sigma_hat, c(0.025, 0.975))
  expect_true(
    sigma_true >= sigma_q[[1]] && sigma_true <= sigma_q[[2]],
    label = paste0("sigma_true = ", sigma_true,
                   ", 95% CI: [", round(sigma_q[[1]], 3),
                   ", ", round(sigma_q[[2]], 3), "]")
  )

  # Extend the existing fit (no new brm() call) to also check that
  # posterior_predict()/posterior_epred() honour resp_trunc()'s bounds --
  # confirms the R-side truncation fix works through the full brms
  # dispatch path, not just the synthetic-prep unit tests below.
  pp <- brms::posterior_predict(fit)
  expect_true(all(sweep(pp, 2, dat$neg_bound, `>=`)),
              label = "posterior_predict draws below their row's trunc(lb=) bound")

  # brms::posterior_epred() itself cannot be used here -- see the matching
  # comment in test-lccdf.R for the confirmed brms-level reason
  # (posterior_epred.brmsprep() unconditionally routes truncated fits to
  # posterior_epred_trunc(), which has no generic handler for ANY custom
  # family and always errors, independent of this package's own fix).
  prep <- brms::prepare_predictions(fit)
  ep   <- posterior_epred_dnorm1(prep)
  expect_true(all(is.finite(ep)))
  expect_true(all(sweep(ep, 2, dat$neg_bound, `>=`)),
              label = "posterior_epred below the row's trunc(lb=) bound")
})

# -----------------------------------------------------------------------
# posterior_predict / posterior_epred truncation correctness (synthetic prep)
# -----------------------------------------------------------------------

test_that("posterior_predict_dnorm1 respects lb (repro of confirmed bug)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 3000, ncol = 1)),
    Y = 0, lb = -3
  )
  set.seed(123)
  draws <- posterior_predict_dnorm1(1, prep)
  expect_true(all(draws >= -3), label = paste0("min draw = ", min(draws)))
})

test_that("posterior_predict_dnorm1 without lb/ub matches untruncated behaviour (fast path)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 5000, ncol = 1)),
    Y = 0
  )
  set.seed(1)
  draws <- posterior_predict_dnorm1(1, prep)
  expect_true(is.numeric(draws) && length(draws) == 5000)
  expect_true(min(draws) < -5)
})

test_that("posterior_predict_dnorm1 draws match the truncated PMF distributionally", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 20000, ncol = 1)),
    Y = 0, lb = -5, ub = 8
  )
  set.seed(42)
  draws <- posterior_predict_dnorm1(1, prep)
  expect_true(all(draws >= -5 & draws <= 8))

  support <- -5:8
  probs <- exp(dnorm1_lpmf_r(support, 4))
  probs <- probs / sum(probs)
  emp <- as.numeric(table(factor(draws, levels = support))) / length(draws)
  expect_lt(max(abs(emp - probs)), 0.02)
})

test_that("posterior_epred_dnorm1 differs from untruncated mean (0) when lb is tight", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(4, nrow = 5, ncol = 1)),
    Y = 0, lb = -3
  )
  epred <- posterior_epred_dnorm1(prep)
  expect_equal(dim(epred), c(5, 1))
  expect_true(all(epred[, 1] > 0), label = paste0("epred = ", epred[1, 1]))
})

test_that("posterior_epred_dnorm1 matches brute-force truncated mean", {
  sigma_val <- 4; lb_val <- -3
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(sigma_val, nrow = 1, ncol = 1)),
    Y = 0, lb = lb_val
  )
  epred <- posterior_epred_dnorm1(prep)

  support <- lb_val:round(10 * sigma_val)
  probs <- exp(dnorm1_lpmf_r(support, sigma_val))
  probs <- probs / sum(probs)
  brute_force_mean <- sum(support * probs)

  expect_equal(epred[1, 1], brute_force_mean, tolerance = 1e-6)
})

test_that("posterior_epred_dnorm1 leaves untruncated observations exactly at 0 (no regression)", {
  prep <- make_synthetic_prep(
    dpars = list(mu = matrix(c(4, 6), nrow = 2, ncol = 2)),
    Y = c(0, 0)
  )
  epred <- posterior_epred_dnorm1(prep)
  expect_equal(epred, 0 * prep$dpars$mu)
})

# -----------------------------------------------------------------------
# ifelse() length-collapse regression -- see test-dlaplace1.R for the full
# explanation. dnorm1_lpmf_r's test (z >= 0) is built from z alone, so a
# scalar z against vector sigma is the direction that would collapse.
# -----------------------------------------------------------------------

test_that("dnorm1_lpmf_r does not collapse to length 1 for scalar z, vector sigma", {
  z <- 2
  sigma <- c(1, 2, 3, 4, 5)
  out <- dnorm1_lpmf_r(z, sigma)
  expect_length(out, length(sigma))
  expect_equal(out, vapply(sigma, function(s) dnorm1_lpmf_r(z, s), numeric(1)))
})
