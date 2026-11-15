# tests/testthat/test-binegbin.R
#
# The fully paired constructor binegbin(): the likelihood, agreement of Stan
# with the R reference, and recovery from simulated data.
#
# The simulated data use one dispersion for both source-exclusive components,
# so the fitted models tie shapexone and shapextwo through one non-linear
# parameter:
#
#   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1
#
# Recovery of two different dispersions is tested in
# test-binegbin-dispersions.R.

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("R brute-force lpmf normalises to 1 across parameter sets", {
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

test_that("moment identities hold (mean, Var(y1), Var(d), Cov)", {
  # With y1 = N_shared + N10, y2 = N_shared + N01, and NB2(m, phi) variance
  # m + m^2/phi: E[y1] = mu + lone; Var(y1) = (mu+mu^2/ss)+(lone+lone^2/sx);
  # Var(d) = 2(lone+lone^2/sx) when lone = ltwo; Cov = mu + mu^2/ss. The
  # tolerances allow for Monte Carlo error with 3e6 draws.
  mu <- 8; lone <- 3; ltwo <- 3; ss <- 2; sx <- 1.5
  set.seed(1); n <- 3e6
  Ns  <- rnbinom(n, size = ss, mu = mu)
  N10 <- rnbinom(n, size = sx, mu = lone)
  N01 <- rnbinom(n, size = sx, mu = ltwo)
  ye <- Ns + N10; yl <- Ns + N01; d <- ye - yl
  expect_equal(mean(ye), mu + lone,                         tolerance = 0.02)
  expect_equal(var(ye),  (mu + mu^2/ss) + (lone + lone^2/sx), tolerance = 0.05)
  expect_equal(var(d),   2 * (lone + lone^2/sx),              tolerance = 0.05)
  expect_equal(cov(ye, yl), mu + mu^2/ss,                   tolerance = 0.1)
})

test_that("binegbin reduces to bipois as shapes, shapex -> Inf (Poisson limit)", {
  # NB2 tends to Poisson as phi -> Inf, with a difference of order 1/phi, so
  # phi = 1e7 is compared at a tolerance of 1e-3.
  mu <- 4; lone <- 2; ltwo <- 3
  ys <- 0:40
  yg <- expand.grid(y1 = ys, y2 = ys)
  lp_nb   <- binegbin_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, 1e7, 1e7)
  lp_pois <- bipois_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo)
  expect_equal(lp_nb, lp_pois, tolerance = 1e-3)
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

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

test_that("Stan binegbin_lpmf matches R brute-force reference across a grid", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu     = c(0.2, 1, 5, 20),
    lone    = c(0.2, 1, 5),
    ltwo    = c(0.2, 1, 5),
    ss     = c(0.5, 2, 50),
    sx     = c(0.5, 2, 50)
  )

  check_one <- function(mu, lone, ltwo, ss, sx) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(function(r, s) binegbin_lpmf(r, mu, lone, ltwo, ss, sx, sx, s, 1L),
                        yg$y1, yg$y2)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, mu, lone, ltwo, ss, sx)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo, grid$ss, grid$sx)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
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
    yg <- expand.grid(y1 = ys, y2 = ys)
    stan_vals <- mapply(
      function(r, s) binegbin_lpmf(r, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], ec[["ss"]], ec[["sx"]], ec[["sx"]], s, 1L),
      yg$y1, yg$y2)
    r_vals <- binegbin_lpmf_r(yg$y1, yg$y2, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], ec[["ss"]], ec[["sx"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# Parameter recovery from simulated data (brms end-to-end)
# -----------------------------------------------------------------------
#
# A. Recovery, intercepts only. A mismatch between the order of the dpars and
#    the Stan lpmf signature exchanges rates or dispersions without an error,
#    and this test detects such a mismatch. The design has no group-level
#    terms and n = 400.
#
# B. Composition with (1 | group). This test checks that the family samples
#    cleanly with a group-level term, and makes no claim about recovery of the
#    dispersions.
#
# The exclusive dispersion is identified through Var(d) alone, and trades off
# against the exclusive rate lambda along a ridge:
#
#   Var(d) = 2 (lambda + lambda^2 / shapex)  =>  shapex = lambda^2 / (V - lambda)
#
# with V = Var(d)/2. The true values below give a variance twice the mean in
# each component (shapes = 8 with mean 8; shapex = 3 with mean 3), which keeps
# 1/phi distinguishable from zero without placing the posterior on the ridge.

BINEGBIN_TRUTH <- list(
  log_mu_int = log(8),
  lone       = 3,   # both exclusive rates; tied to one value via nlf(~ lamx)
  shapes     = 8,
  shapex     = 3
)

# brms names the intercept of mu b_Intercept and the intercept of shapes
# b_shapes_Intercept. The exclusive dispersion is the non-linear parameter
# `shapexx`, with intercept b_shapexx_Intercept. Values are on the log scale.
BINEGBIN_DRAWS_TRUTH <- c(
  b_Intercept         = BINEGBIN_TRUTH$log_mu_int,
  b_lamx_Intercept    = log(BINEGBIN_TRUTH$lone),
  b_shapes_Intercept  = log(BINEGBIN_TRUTH$shapes),
  b_shapexx_Intercept = log(BINEGBIN_TRUTH$shapex)
)

# Intercepts-only simulation, used by the recovery test and the coverage
# assessment.
binegbin_sim <- function(seed, n = 400L) {
  set.seed(seed)
  n_shared <- rnbinom(n, size = BINEGBIN_TRUTH$shapes,
                      mu = exp(BINEGBIN_TRUTH$log_mu_int))
  n1 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  n2 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  data.frame(y1 = n_shared + n1, y2 = n_shared + n2)
}

# The same true values with a group effect on mu, for the composition test.
BINEGBIN_RE_SD <- 0.3
binegbin_sim_re <- function(seed, n_group = 8L, n_per_group = 25L) {
  set.seed(seed)
  n      <- n_group * n_per_group
  group <- rep(seq_len(n_group), each = n_per_group)
  mu_i   <- exp(BINEGBIN_TRUTH$log_mu_int +
                  BINEGBIN_RE_SD * rnorm(n_group)[group])
  n_shared <- rnbinom(n, size = BINEGBIN_TRUTH$shapes, mu = mu_i)
  n1 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  n2 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  data.frame(y1 = n_shared + n1, y2 = n_shared + n2, group = factor(group))
}

# `mu` requires its own formula. Under nl = TRUE, the right-hand side of the
# main formula is a non-linear expression for mu, so `bf(y1 | vint(y2) ~ 1,
# nl = TRUE)` without `mu ~ 1` fixes the shared rate at exp(1) without an
# error.
binegbin_fit <- function(dat, mu_re = FALSE, ...) {
  f <- if (mu_re) {
    brms::bf(
      y1 | vint(y2) ~ 1,
      mu ~ 1 + (1 | group),
      brms::nlf(lambdaone ~ lamx),
      brms::nlf(lambdatwo ~ lamx),
      brms::nlf(shapexone ~ shapexx),
      brms::nlf(shapextwo ~ shapexx),
      lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE
    )
  } else {
    brms::bf(
      y1 | vint(y2) ~ 1,
      mu ~ 1,
      brms::nlf(lambdaone ~ lamx),
      brms::nlf(lambdatwo ~ lamx),
      brms::nlf(shapexone ~ shapexx),
      brms::nlf(shapextwo ~ shapexx),
      lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE
    )
  }
  # The brms defaults give `lamx`, `shapes`, and `shapexx` flat improper
  # priors. Without the priors below, this fit gave a divergent transition.
  # The priors are centred at zero and the true values are positive (log 8 =
  # 2.08, log 3 = 1.10), so each prior pulls away from the true value.
  # normal(0, 2) on a log dispersion spans 0.02 to 50 at 95%, and
  # normal(0, 3) on a log rate spans 0.003 to 357.
  prm <- c(
    brms::prior(normal(0, 3), class = "Intercept"),
    brms::prior(normal(0, 3), class = "b", nlpar = "lamx"),
    brms::prior(normal(0, 2), class = "Intercept", dpar = "shapes"),
    brms::prior(normal(0, 2), class = "b", nlpar = "shapexx")
  )

  suppressMessages(brms::brm(
    f,
    family   = binegbin(),
    stanvars = binegbin_stanvars(),
    data     = dat,
    prior    = prm,
    backend  = "rstan",
    chains   = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 20260705,
    refresh  = 0,
    control  = list(adapt_delta = 0.99),
    ...
  ))
}

# A. Recovery -------------------------------------------------------------

test_that("binegbin recovers all five dpars (intercepts only)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  dat   <- binegbin_sim(20260705)
  fit   <- binegbin_fit(dat)
  draws <- as.data.frame(fit)

  # Smoke gates at a wide interval; see helper-coverage.R. Calibration is
  # assessed by the coverage test below.
  for (p in names(BINEGBIN_DRAWS_TRUTH)) {
    expect_true(
      recovery_ok(draws, BINEGBIN_DRAWS_TRUTH[[p]], p),
      label = paste0("recovery of ", p)
    )
  }

  # The observed Var(d) must lie inside the 99% posterior predictive interval.
  # Var(d) depends on the exclusive dispersion. Both counts are drawn from the
  # joint model, because subtracting the recorded y2 from marginal draws of y1
  # would overstate Var(d) by 2 Var(N_shared).
  pairs <- posterior_predict_y1_and_y2(fit, ndraws = 200)
  drep  <- pairs[, , "y1"] - pairs[, , "y2"]
  q <- quantile(apply(drep, 1, var), c(0.005, 0.995), names = FALSE)
  var_d_obs <- var(dat$y1 - dat$y2)
  expect_true(
    var_d_obs >= q[[1]] && var_d_obs <= q[[2]],
    label = paste0("observed Var(d) = ", round(var_d_obs, 2),
                   " outside predictive 99% [", round(q[[1]], 2), ", ",
                   round(q[[2]], 2), "]")
  )

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# B. Composition with a group-level term --------------------------------

test_that("binegbin composes with a group-level term on mu and samples cleanly", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  fit   <- binegbin_fit(binegbin_sim_re(20260706), mu_re = TRUE)
  draws <- as.data.frame(fit)

  # With 8 groups, the group-level SD is weakly identified and trades off
  # against the dispersions, so the test makes no recovery assertion on the
  # dispersions.
  expect_true("sd_group__Intercept" %in% names(draws))
  sd_draws <- draws[["sd_group__Intercept"]]
  expect_true(all(is.finite(sd_draws)))
  expect_true(all(sd_draws > 0))

  # A broken group-level term would give a median SD outside (0.02, 3).
  expect_true(
    median(sd_draws) > 0.02 && median(sd_draws) < 3,
    label = paste0("median sd_group = ", round(median(sd_draws), 3))
  )

  # Every dpar intercept is present and finite.
  for (p in names(BINEGBIN_DRAWS_TRUTH)) {
    expect_true(p %in% names(draws), label = paste0(p, " present"))
    expect_true(all(is.finite(draws[[p]])), label = paste0(p, " finite"))
  }

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.02, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# Calibration --------------------------------------------------------------

test_that("binegbin posterior intervals are calibrated (coverage over replicates)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()
  skip_unless_coverage()

  R           <- 10L
  level       <- 0.9
  floor_count <- coverage_floor(R, level)   # 6 of 10 at alpha = 0.01

  # The intercepts-only model is compiled once, and update() re-runs sampling.
  # Replicate seeds differ from the seed of the recovery fit.
  fit0 <- binegbin_fit(binegbin_sim(20260705))

  # The replicates run 2 chains of 1500 post-warmup draws, enough for the 5%
  # and 95% quantiles.
  cov <- coverage_recovery(
    fit0,
    sim    = function(i) binegbin_sim(20260800 + i),
    truths = BINEGBIN_DRAWS_TRUTH,
    R      = R,
    level  = level,
    chains = 2, iter = 2000, warmup = 500
  )

  message("coverage at nominal ", level, " over ", R, " replicates:")
  message(paste0("  ", names(cov), ": ", cov * R, "/", R, collapse = "\n"))

  for (p in names(cov)) {
    expect_gte(
      cov[[p]] * R, floor_count,
      label = paste0("coverage of ", p, " (", cov[[p]] * R, "/", R,
                     "; floor ", floor_count, "/", R,
                     " at nominal ", level, ")")
    )
  }
})
