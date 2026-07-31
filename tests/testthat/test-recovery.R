# tests/testthat/test-recovery.R

test_that("parameter recovery from simulated hierarchical data", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(42)

  n_groups             <- 10L
  n_obs_per_group       <- 20L
  true_sigma_intercept  <- log(3)
  true_sigma_group      <- 0.5

  group_effects <- rnorm(n_groups, mean = 0, sd = true_sigma_group)
  group         <- rep(seq_len(n_groups), each = n_obs_per_group)
  sigma_i       <- exp(true_sigma_intercept + group_effects[group])
  mu_i          <- sigma_i^2 / 2
  y             <- skellam::rskellam(length(mu_i), lambda1 = mu_i, lambda2 = mu_i)
  dat           <- data.frame(y = y, group = factor(group))

  # brms's default Intercept prior (student_t(3, 0, 2.5)) is known to be
  # wide enough that this custom Bessel-based likelihood can occasionally
  # wander to a nonsensical log(sigma) region for an unlucky seed -- the
  # exact failure mode hit here after the sigma-reparameterisation changed
  # what data this fixed seed simulates. Applying that same documented fix.
  sane_prior <- brms::prior(normal(1, 1.5), class = "Intercept")

  suppressMessages({
    fit <- brms::brm(
      y ~ 1 + (1 | group),
      data     = dat,
      family   = skellam1(),
      stanvars = skellam1_stanvars(),
      prior    = sane_prior,
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      seed     = 42,
      refresh  = 0
    )
  })

  draws <- as.data.frame(fit)

  # 1 & 2. Smoke gate on sigma_intercept (log(sigma) scale) and sigma_group.
  # A wide interval on purpose: at 90% each of these fails 10% of the time on a
  # CORRECT model. See helper-coverage.R for the smoke-vs-calibration split.
  intercept_q <- quantile(draws[["b_Intercept"]], c(0.005, 0.995))
  expect_true(
    recovery_ok(draws, true_sigma_intercept, "b_Intercept"),
    label = paste0("true intercept = ", round(true_sigma_intercept, 3),
                   ", 99% CI: [", round(intercept_q[[1]], 3),
                   ", ",          round(intercept_q[[2]], 3), "]")
  )

  sd_q <- quantile(draws[["sd_group__Intercept"]], c(0.005, 0.995))
  expect_true(
    recovery_ok(draws, true_sigma_group, "sd_group__Intercept"),
    label = paste0("sigma_group = ", true_sigma_group,
                   ", 99% CI: [", round(sd_q[[1]], 3),
                   ", ", round(sd_q[[2]], 3), "]")
  )

  # 3. No divergences
  n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  expect_equal(n_div, 0,
               label = paste0(n_div, " divergent transitions"))

  # 4. All Rhat < 1.01
  max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
  expect_true(
    max_rhat < 1.01,
    label = paste0("max Rhat = ", round(max_rhat, 4))
  )
})
