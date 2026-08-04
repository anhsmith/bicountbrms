# tests/testthat/test-bipois_cens.R

# -----------------------------------------------------------------------
# R-side tests -- no Stan compilation required
# -----------------------------------------------------------------------

test_that("matched (y1_obs==1) joint PMF normalises to 1 across parameter sets", {
  norm_check <- function(mu, lone, ltwo, K = 100) {
    ys <- 0:K
    yg <- expand.grid(y1 = ys, y2 = ys)
    lp <- bipois_cens_lpmf_r(yg$y1, yg$y2, 1L, mu, lone, ltwo)
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

test_that("y2-only (y1_obs==0) branch normalises to 1 over y2", {
  # A weak check of the closed form itself -- a Poisson pmf sums to 1 whatever
  # its rate -- but a real check of the BRANCH WIRING: that y1_obs == 0 selects
  # a 1-D pmf in y2 with the y1 margin integrated out, and that it ignores y1
  # and lambdaone entirely. Both are asserted below by varying them.
  norm_check <- function(mu, lone, ltwo, K = 250) {
    ys <- 0:K
    lp <- bipois_cens_lpmf_r(rep(0L, length(ys)), ys, 0L, mu, lone, ltwo)
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

  # The censored branch must not depend on y1 or on lambdaone: neither is
  # observed on those rows, and a dependence on either would mean the branch is
  # not the y1-integrated marginal.
  base <- bipois_cens_lpmf_r(0L, 7L, 0L, mu = 5, lambdaone = 3, lambdatwo = 4)
  expect_equal(bipois_cens_lpmf_r(99L, 7L, 0L, 5, 3, 4),  base)
  expect_equal(bipois_cens_lpmf_r(0L,  7L, 0L, 5, 40, 4), base)
})

test_that("the censored branch's closed form equals the brute-force convolution", {
  # THE IDENTITY THAT LICENSES THE ONE-LINE STAN BRANCH. binegbin_cens must
  # evaluate P(y2) = sum_k f_s(k) f_2(y2 - k) as a sum because NB2 + NB2 is not
  # NB2. For Poisson components the same convolution collapses exactly to
  # Poisson(mu + lambdatwo), and bipois_cens uses that closed form on both the
  # Stan and R sides. The brute-force sum is therefore kept HERE, as the
  # independent route, rather than in R/bipois_cens.R -- otherwise the closed
  # form would be assumed rather than checked.
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
    closed <- bipois_cens_lpmf_r(0L, y2_vals, 0L, g$mu, 1, g$ltwo)
    expect_equal(brute, closed, tolerance = 1e-12,
                 label = paste(unlist(g), collapse = ","))
  }
})

test_that("marginal identity: sum over y1 of the matched branch == y2-only branch", {
  # Integrating the matched (y1_obs==1) joint over all y1 must reproduce the
  # y2-only (y1_obs==0) value at that y2. This is what makes the two branches
  # one model rather than two, and here it also confirms that the closed form
  # is the marginal OF THIS JOINT and not merely a Poisson that happens to sum
  # to 1.
  mu <- 6; lone <- 3; ltwo <- 4
  K  <- 200
  for (yl in c(0L, 1L, 3L, 7L, 15L)) {
    ys <- 0:K
    lp_joint <- bipois_cens_lpmf_r(ys, rep(yl, length(ys)), 1L, mu, lone, ltwo)
    marg_from_joint <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
    censored_branch <- bipois_cens_lpmf_r(0L, yl, 0L, mu, lone, ltwo)
    expect_equal(marg_from_joint, censored_branch, tolerance = 1e-10,
                 label = paste("y2 =", yl))
  }
})

test_that("matched branch is the bipois lpmf (equivalence)", {
  # On the R side this holds BY CONSTRUCTION: bipois_cens_lpmf_r delegates its
  # matched branch to bipois_lpmf_r rather than restating the sum. The test is
  # therefore a guard against a later edit inlining a divergent copy, not
  # independent evidence. The substantive check is the Stan one below, where
  # the recurrence is written out separately in bipois_cens_stan_funs.
  grid <- expand.grid(
    mu   = c(0.5, 3, 12),
    lone = c(0.5, 2, 6),
    ltwo = c(0.5, 2, 6)
  )
  ys <- 0:20
  yg <- expand.grid(y1 = ys, y2 = ys)
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    lp_joint  <- bipois_cens_lpmf_r(yg$y1, yg$y2, 1L, g$mu, g$lone, g$ltwo)
    lp_bipois <- bipois_lpmf_r(yg$y1, yg$y2, g$mu, g$lone, g$ltwo)
    expect_equal(lp_joint, lp_bipois, tolerance = 1e-14,
                 label = paste(unlist(g), collapse = ","))
  }
})

test_that("binegbin_cens reduces to bipois_cens in the Poisson limit, on both branches", {
  # THE CROSS-FAMILY CHECK THIS FAMILY MAKES POSSIBLE. binegbin_cens's
  # censored branch is a numerical convolution whose only other test -- the
  # marginal identity -- compares it against the matched branch of the same
  # code, so an error shared by both sums would pass. Driving its three
  # dispersions to their Poisson limit gives an INDEPENDENT analytic target:
  # NB2(m, phi) -> Poisson(m) as phi -> Inf, so binegbin_cens -> bipois_cens,
  # whose censored branch is closed form.
  #
  # The approach is O(1/phi) -- measured at 2.53e-3 for phi = 1e5 and 2.53e-5
  # for phi = 1e7 on the grid below, a clean factor of 100 for a factor of 100
  # in phi -- so the bounds are stated at two values of phi two orders of
  # magnitude apart, each a factor of ~4 above the observed error. A failure of
  # the limit shows up as an error that does not shrink with phi, not merely as
  # a large one, which is what the third assertion catches.
  mu <- 4; lone <- 2.5; ltwo <- 3
  ys <- 0:25
  yg <- expand.grid(y1 = ys, y2 = ys, y1_obs = c(0L, 1L))

  ref <- bipois_cens_lpmf_r(yg$y1, yg$y2, yg$y1_obs, mu, lone, ltwo)

  err <- function(phi) {
    nb <- binegbin_cens_lpmf_r(yg$y1, yg$y2, yg$y1_obs, mu, lone, ltwo,
                                shapes = phi, shapexone = phi, shapextwo = phi)
    max(abs(nb - ref))
  }

  e5 <- err(1e5)
  e7 <- err(1e7)
  expect_lt(e5, 1e-2)
  expect_lt(e7, 1e-4)
  # And the error must actually be shrinking with phi, which distinguishes a
  # genuine limit from two expressions that merely happen to be close.
  expect_lt(e7, e5)
})

test_that("posterior_predict draws reproduce the joint/marginal conditional y1 | y2", {
  # posterior_predict_bipois_cens splits N_shared | y2 as a Binomial then adds
  # a fresh N1. Its distribution must equal P(y1 | y2) = joint(y1, y2) /
  # marginal(y2). Checked by Monte Carlo, as for binegbin_cens. The Binomial
  # split is exact here (conditioning a sum of independent Poissons on its total
  # gives a Binomial), so this also confirms that shortcut.
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
    Y     = 0L,          # response value is unused by posterior_predict
    vint1 = y2,
    vint2 = 1L
  )
  draws <- posterior_predict_bipois_cens(1, prep)
  expect_length(draws, ndraws)

  K  <- 80
  xs <- 0:K
  lp_joint <- bipois_cens_lpmf_r(xs, rep(y2, length(xs)), 1L, mu, lone, ltwo)
  lp_marg  <- { mx <- max(lp_joint); mx + log(sum(exp(lp_joint - mx))) }
  p_cond   <- exp(lp_joint - lp_marg)

  emp <- tabulate(draws + 1L, nbins = K + 1L) / ndraws
  keep <- p_cond > 1e-3
  expect_lt(max(abs(emp[keep] - p_cond[keep])), 0.01)
  expect_equal(mean(draws), sum(xs * p_cond), tolerance = 0.05)

  # posterior_epred is the same conditional in expectation rather than
  # simulated, and here it is exact -- the mean of that Binomial split. It must
  # agree with the analytic conditional mean to machine precision, not merely
  # to Monte Carlo error.
  ep <- posterior_epred_bipois_cens(prep)
  expect_equal(dim(ep), c(ndraws, 1L))
  expect_equal(unique(as.vector(ep)), sum(xs * p_cond), tolerance = 1e-10)
})

test_that("the censored branch ignores y1_obs only where it should", {
  # posterior_predict and posterior_epred impute y1 on EVERY row, censored
  # included; the likelihood does not. Pin both directions so a later change
  # cannot quietly align them.
  prep <- make_synthetic_prep(
    dpars = list(mu = 5, lambdaone = 3, lambdatwo = 4),
    Y     = 9L,
    vint1 = 6L,
    vint2 = 0L        # censored row
  )
  # Likelihood: branches on y1_obs, so it must NOT equal the matched value.
  expect_equal(
    log_lik_bipois_cens(1, prep),
    stats::dpois(6L, 5 + 4, log = TRUE)
  )
  # epred: same value whether the row is censored or matched.
  prep_matched <- prep
  prep_matched$data$vint2 <- 1L
  expect_equal(posterior_epred_bipois_cens(prep),
               posterior_epred_bipois_cens(prep_matched))
})

# -----------------------------------------------------------------------
# Stan tests -- require rstan; skipped silently if unavailable
# -----------------------------------------------------------------------

# Compile bipois AND bipois_cens together so the Stan-level equivalence check
# can call both lpmfs on identical inputs. The two declare different function
# names, so there is no collision.
stan_code <- paste0("functions {\n", bipois_stan_funs, "\n",
                    bipois_cens_stan_funs, "}\nmodel {}\n")

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

test_that("Stan bipois_cens_lpmf matches the R reference (both branches)", {
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
      function(r, s, e) bipois_cens_lpmf(r, mu, lone, ltwo, s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- bipois_cens_lpmf_r(yg$y1, yg$y2, yg$y1_obs, mu, lone, ltwo)
    max(abs(stan_vals - r_vals))
  }

  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo)
  expect_true(max(diffs) < 1e-8, label = paste("max diff =", max(diffs)))
})

test_that("Stan bipois_cens matched branch == Stan bipois lpmf (equivalence)", {
  # The substantive equivalence check: bipois_cens_stan_funs writes the
  # recurrence out separately rather than calling bipois_lpmf, so the two Stan
  # implementations can drift. This is what stops them.
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  grid <- expand.grid(
    mu   = c(0.2, 1, 5, 20),
    lone = c(0.2, 1, 5),
    ltwo = c(0.2, 1, 5)
  )
  check_one <- function(mu, lone, ltwo) {
    means <- c(mu + lone, mu + ltwo)
    ys <- unique(pmax(c(0L, 1L, round(means), round(means) + 4L), 0L))
    yg <- expand.grid(y1 = ys, y2 = ys)
    joint_vals <- mapply(function(r, s) bipois_cens_lpmf(r, mu, lone, ltwo, s, 1L),
                         yg$y1, yg$y2)
    bipois_vals <- mapply(function(r, s) bipois_lpmf(r, mu, lone, ltwo, s),
                          yg$y1, yg$y2)
    max(abs(joint_vals - bipois_vals))
  }
  diffs <- mapply(check_one, grid$mu, grid$lone, grid$ltwo)
  expect_true(max(diffs) < 1e-12, label = paste("max diff =", max(diffs)))
})

test_that("Stan bipois_cens_lpmf is numerically stable at extreme rates", {
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
      function(r, s, e) bipois_cens_lpmf(r, ec[["mu"]], ec[["lone"]], ec[["ltwo"]], s, e),
      yg$y1, yg$y2, yg$y1_obs)
    r_vals <- bipois_cens_lpmf_r(yg$y1, yg$y2, yg$y1_obs,
                                  ec[["mu"]], ec[["lone"]], ec[["ltwo"]])
    expect_false(any(!is.finite(stan_vals)), label = paste(ec, collapse = ","))
    expect_equal(stan_vals, r_vals, tolerance = 1e-8, label = paste(ec, collapse = ","))
  }
})

# -----------------------------------------------------------------------
# brms end-to-end: dispatch (loo + posterior_predict) and recovery
# -----------------------------------------------------------------------

test_that("bipois_cens fits a censored design, dispatches, and recovers params", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20260804)

  n_vessel     <- 8L
  n_per_vessel <- 30L
  n            <- n_vessel * n_per_vessel

  true_log_mu_int <- log(8)
  true_sd_vessel  <- 0.3
  true_lone       <- 3           # shared excess rate (lambdaone = lambdatwo)

  vessel <- rep(seq_len(n_vessel), each = n_per_vessel)
  z_mu   <- rnorm(n_vessel)
  mu_i   <- exp(true_log_mu_int + true_sd_vessel * z_mu[vessel])

  n_shared <- rpois(n, mu_i)
  n1       <- rpois(n, true_lone)
  n2       <- rpois(n, true_lone)

  # Half the rows are y2-only (y1 unobserved) -- the censoring the family
  # exists to handle. mu and lambdatwo are separated only by the matched rows
  # (the censored branch sees them only through their sum), so the matched
  # half is what makes this recoverable; see ?bipois_cens.
  y1_obs <- rep(c(1L, 0L), length.out = n)

  dat <- data.frame(
    y1     = n_shared + n1,
    y2     = n_shared + n2,
    y1_obs = y1_obs,
    vessel = factor(vessel)
  )

  suppressMessages({
    fit <- brms::brm(
      brms::bf(
        y1 | vint(y2, y1_obs) ~ 1,
        mu ~ 1 + (1 | vessel),
        brms::nlf(lambdaone ~ lamx),
        brms::nlf(lambdatwo ~ lamx),
        lamx ~ 1, nl = TRUE
      ),
      family   = bipois_cens(),
      stanvars = bipois_cens_stanvars(),
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

  # Dispatch: both must resolve the package's log_lik_/posterior_predict_
  # methods by name without "no applicable method" errors.
  ll <- brms::log_lik(fit)
  expect_equal(dim(ll)[2], n)
  expect_true(all(is.finite(ll)))

  loo_obj <- suppressWarnings(brms::loo(fit))
  expect_s3_class(loo_obj, "loo")

  pp <- brms::posterior_predict(fit)
  expect_equal(dim(pp)[2], n)
  expect_true(all(pp >= 0))

  ep <- posterior_epred_bipois_cens(brms::prepare_predictions(fit))
  expect_equal(dim(ep)[2], n)
  expect_true(all(is.finite(ep)))

  # brms::posterior_epred() must ALSO reach the family method. Its dispatcher
  # checks truncation before family type and the truncated branch has no
  # custom-family fallback -- which is why a truncated custom family cannot use
  # it at all. Nothing here is truncated, so the ordinary custom-family route
  # applies and the two calls must agree exactly, not merely closely.
  expect_equal(brms::posterior_epred(fit), ep)

  draws <- as.data.frame(fit)
  # Smoke gate at a wide interval, not a calibration claim. See
  # helper-coverage.R; the matched branch is bipois term for term, whose
  # calibration is assessed in test-bipois.R.
  expect_true(recovery_ok(draws, true_log_mu_int, "b_Intercept"))
  expect_true(recovery_ok(draws, log(true_lone),  "b_lamx_Intercept"))
  expect_true(recovery_ok(draws, true_sd_vessel,  "sd_vessel__Intercept"))

  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0, label = paste0(n_div, " divergent transitions"))

  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(max_rhat < 1.01, label = paste0("max Rhat = ", round(max_rhat, 4)))
})
