# tests/testthat/test-binegbin.R

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
  # For y1 = N_shared + N10, y2 = N_shared + N01, NB2(m, phi) with
  # var = m + m^2/phi: E[y1] = mu + lone; Var(y1) = (mu+mu^2/ss)+(lone+lone^2/sx);
  # Var(d) = 2(lone+lone^2/sx) when lone=ltwo, sx shared; Cov = mu + mu^2/ss.
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
  # NB2 -> Poisson as phi -> Inf with O(1/phi) residual; a large-but-finite phi
  # leaves a small tail difference, so use a generous phi and a modest tolerance.
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
if (requireNamespace("rstan", quietly = TRUE)) {
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
    stan_vals <- mapply(function(r, s) binegbin_lpmf(r, mu, lone, ltwo, ss, sx, s),
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
      function(r, s) binegbin_lpmf(r, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], ec[["ss"]], ec[["sx"]], s),
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
# Two tests, deliberately doing different jobs.
#
# A. RECOVERY (intercepts only). Validates the family's parameterisation and,
#    critically, the POSITIONAL correspondence between the dpars vector and the
#    Stan lpmf signature -- a mismatch there does not error, it silently swaps
#    which rate or dispersion governs which component. Kept well conditioned on
#    purpose: no group-level terms, n = 400, moderate overdispersion.
#
# B. COMPOSITION (adds (1 | vessel)). Checks only that the family composes with
#    brms's group-level machinery and samples cleanly. It makes NO recovery
#    claim about the dispersions -- see below.
#
# WHY THE SPLIT. One test used to do both jobs at once, with a vessel effect on
# 8 levels and severe overdispersion (the shared component's variance was 40
# against a mean of 8, i.e. 80% of it overdispersion). Three variance channels
# -- the group-level SD, the shared dispersion, the private dispersion -- then
# competed for the same residual, and shapex lost. shapex is identified only
# through the difference's variance, and with lambda also unknown the two trade
# off along a ridge:
#
#   Var(d) = 2 (lambda + lambda^2 / shapex)  =>  shapex = lambda^2 / (V - lambda)
#
# with V = Var(d)/2. At fixed V, every lambda implies a shapex. On the old
# seeded dataset Var(d) came out 20.7 against the 18 the truth implied -- an
# ordinary ~1.5 SD draw at n = 200 -- and the posterior slid along the ridge to
# (lambda, shapex) = (2.49, 0.795) against a truth of (3, 1.5). A grid MLE of
# the exact likelihood agreed to two decimals, so that was the likelihood
# tracking the data, not a fault: the test was asserting a property of one
# DRAW, not of the estimator. Calibration of the estimator is what the coverage
# test at the bottom of this file measures.
#
# The truths below put ~50% of each component's variance into overdispersion
# (shapes = 8 gives var 16 vs mean 8; shapex = 3 gives var 6 vs mean 3). That is
# identifiable from both directions -- far enough from the Poisson limit that
# 1/phi is clearly non-zero, far enough from the extreme that the ridge does not
# dominate.

BINEGBIN_TRUTH <- list(
  log_mu_int = log(8),
  lone       = 3,   # both private rates; tied to one value via nlf(~ lamx)
  shapes     = 8,
  shapex     = 3
)

# shapes/shapex are log-linked, so brms reports them as b_<dpar>_Intercept.
# "mu" is the family's canonical dpar, so brms drops its infix (b_Intercept).
BINEGBIN_DRAWS_TRUTH <- c(
  b_Intercept        = BINEGBIN_TRUTH$log_mu_int,
  b_lamx_Intercept   = log(BINEGBIN_TRUTH$lone),
  b_shapes_Intercept = log(BINEGBIN_TRUTH$shapes),
  b_shapex_Intercept = log(BINEGBIN_TRUTH$shapex)
)

# Intercepts-only generative draw, shared by the recovery test and the coverage
# assessment so the two cannot drift apart.
binegbin_sim <- function(seed, n = 400L) {
  set.seed(seed)
  n_shared <- rnbinom(n, size = BINEGBIN_TRUTH$shapes,
                      mu = exp(BINEGBIN_TRUTH$log_mu_int))
  n1 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  n2 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  data.frame(y1 = n_shared + n1, y2 = n_shared + n2)
}

# Same truth plus a vessel effect on mu, for the composition test only.
BINEGBIN_RE_SD <- 0.3
binegbin_sim_re <- function(seed, n_vessel = 8L, n_per_vessel = 25L) {
  set.seed(seed)
  n      <- n_vessel * n_per_vessel
  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  mu_i   <- exp(BINEGBIN_TRUTH$log_mu_int +
                  BINEGBIN_RE_SD * rnorm(n_vessel)[vessel])
  n_shared <- rnbinom(n, size = BINEGBIN_TRUTH$shapes, mu = mu_i)
  n1 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  n2 <- rnbinom(n, size = BINEGBIN_TRUTH$shapex, mu = BINEGBIN_TRUTH$lone)
  data.frame(y1 = n_shared + n1, y2 = n_shared + n2, vessel = factor(vessel))
}

# NOTE: `mu` MUST get its own explicit formula. Under nl = TRUE the main
# formula's right-hand side is a NON-LINEAR EXPRESSION for mu, not a request for
# an intercept -- so `bf(y1 | vint(y2) ~ 1, nl = TRUE)` without a `mu ~ 1` term
# generates `mu[n] = exp(1);`, pinning the shared rate at e and never estimating
# it. It does not error; brms just fits a model you did not ask for, and with the
# level forced into the excess rates the sampler crawls (observed: 150+ minutes
# without finishing one fit). With `mu ~ 1` present the generated code is
# `mu += Intercept; mu = exp(mu);` as intended.
binegbin_fit <- function(dat, mu_re = FALSE, ...) {
  f <- if (mu_re) {
    brms::bf(
      y1 | vint(y2) ~ 1,
      mu ~ 1 + (1 | vessel),
      brms::nlf(lambdaone ~ lamx),
      brms::nlf(lambdatwo ~ lamx),
      lamx ~ 1, shapes ~ 1, shapex ~ 1, nl = TRUE
    )
  } else {
    brms::bf(
      y1 | vint(y2) ~ 1,
      mu ~ 1,
      brms::nlf(lambdaone ~ lamx),
      brms::nlf(lambdatwo ~ lamx),
      lamx ~ 1, shapes ~ 1, shapex ~ 1, nl = TRUE
    )
  }
  # PRIORS MATTER HERE. get_prior() on this model shows brms's defaults leave
  # `lamx`, `shapes` and `shapex` with FLAT IMPROPER priors -- only mu gets a
  # student_t. Those three are the weakly-identified ones, so unregularised they
  # let the chains wander into the flat tail: without these priors this fit
  # produced a divergent transition and max Rhat 1.0101 (against a 1.01 gate),
  # while every recovery assertion still passed. The package documentation tells
  # users to set priors, so the test should exercise that configuration rather
  # than one nobody should use.
  #
  # All four are centred at ZERO while the truths are positive (log 8 = 2.08,
  # log 3 = 1.10), so every prior pulls AWAY from the truth. Recovery therefore
  # remains a genuine test of what the data carries, not a prior echo. They are
  # weak: normal(0, 2) on a log dispersion spans phi in [0.02, 50] at 95%, and
  # normal(0, 3) on a log rate spans [0.003, 357].
  prm <- c(
    brms::prior(normal(0, 3), class = "Intercept"),
    brms::prior(normal(0, 3), class = "b", nlpar = "lamx"),
    brms::prior(normal(0, 2), class = "Intercept", dpar = "shapes"),
    brms::prior(normal(0, 2), class = "Intercept", dpar = "shapex")
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

# A. RECOVERY -----------------------------------------------------------

test_that("binegbin recovers all five dpars (intercepts only)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  dat   <- binegbin_sim(20260705)
  fit   <- binegbin_fit(dat)
  draws <- as.data.frame(fit)

  # A wide interval, so a correctly-specified model does not fail by chance --
  # this is a gate against gross error, not a calibration claim. Calibration is
  # asserted by the coverage test below. See helper-coverage.R.
  for (p in names(BINEGBIN_DRAWS_TRUTH)) {
    expect_true(
      recovery_ok(draws, BINEGBIN_DRAWS_TRUTH[[p]], p),
      label = paste0("recovery of ", p)
    )
  }

  # Independent of the point checks: does the fitted model reproduce the spread
  # shapex actually governs? Var(d) is what identifies it, so a mis-wired
  # dispersion shows up here even if a point interval happened to cover.
  yrep  <- brms::posterior_predict(fit, ndraws = 200)
  y2mat <- matrix(dat$y2, nrow = nrow(yrep), ncol = ncol(yrep), byrow = TRUE)
  q <- quantile(apply(yrep - y2mat, 1, var), c(0.005, 0.995), names = FALSE)
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
  expect_true(max_rhat < 1.01, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# B. COMPOSITION WITH A GROUP-LEVEL TERM --------------------------------

test_that("binegbin composes with a group-level term on mu and samples cleanly", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  fit   <- binegbin_fit(binegbin_sim_re(20260706), mu_re = TRUE)
  draws <- as.data.frame(fit)

  # Deliberately NO dispersion-recovery assertions. With 8 groups the
  # group-level SD is barely identified, and its slack is absorbed by
  # shapes/shapex; asserting tight recovery here tests the draw, not the code.
  # The random-effects machinery is brms's own -- what needs checking is that
  # this custom family composes with it and samples cleanly.
  expect_true("sd_vessel__Intercept" %in% names(draws))
  sd_draws <- draws[["sd_vessel__Intercept"]]
  expect_true(all(is.finite(sd_draws)))
  expect_true(all(sd_draws > 0))

  # A generous bound: a broken group-level term would blow up or collapse, and
  # either would fall outside this.
  expect_true(
    median(sd_draws) > 0.02 && median(sd_draws) < 3,
    label = paste0("median sd_vessel = ", round(median(sd_draws), 3))
  )

  # The dpars must still all be present and finite under the added structure.
  for (p in names(BINEGBIN_DRAWS_TRUTH)) {
    expect_true(p %in% names(draws), label = paste0(p, " present"))
    expect_true(all(is.finite(draws[[p]])), label = paste0(p, " finite"))
  }

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.01, label = paste0("max Rhat = ", round(max_rhat, 4)))
})

# CALIBRATION ------------------------------------------------------------

test_that("binegbin posterior intervals are calibrated (coverage over replicates)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()
  skip_unless_coverage()

  R           <- 10L
  level       <- 0.9
  floor_count <- coverage_floor(R, level)   # 6 of 10 at alpha = 0.01

  # Compile once on the well-conditioned intercepts-only model; update() re-runs
  # sampling only. Replicate seeds are offset from the recovery fit's so the
  # first replicate is not that same dataset.
  fit0 <- binegbin_fit(binegbin_sim(20260705))

  # Shorter chains for the replicates: 3000 draws is ample for a 5%/95%
  # quantile, and cuts each refit to a fraction of the parent fit's cost.
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
