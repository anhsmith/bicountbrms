# tests/testthat/test-y1obs-indicator.R
#
# The check on the observation indicator y1_obs. brms binds vint() terms by
# position, so vint(y1_obs, y2) would pass the second count as the indicator.
# The Stan lpmf of each family rejects an indicator other than 0 or 1, and
# the R references stop with an error that names the correct order.

# ---------------------------------------------------------------------------
# Fast tests -- no Stan compilation
# ---------------------------------------------------------------------------

test_that("each Stan function string rejects an indicator other than 0 or 1", {
  for (scode in list(bipois_stan_funs, binegbin_stan_funs)) {
    expect_true(grepl("if (y1_obs != 0 && y1_obs != 1)", scode, fixed = TRUE))
    expect_true(grepl("reject(\"y1_obs must be 0 or 1; found \", y1_obs,",
                      scode, fixed = TRUE))
    # The check precedes the branch on y1_obs.
    expect_lt(regexpr("reject(", scode, fixed = TRUE),
              regexpr("if (y1_obs == 1)", scode, fixed = TRUE))
  }
})

test_that("the R references stop on an indicator other than 0 or 1", {
  msg <- "y1_obs must be 0 or 1; found 2. Check the order: vint(y2, y1_obs)."
  expect_error(bipois_lpmf_r(3L, 4L, 5, 3, 4, y1_obs = 2L),
               msg, fixed = TRUE)
  expect_error(binegbin_lpmf_r(3L, 4L, 5, 3, 4, 2, 1.5, y1_obs = 2L),
               msg, fixed = TRUE)

  # A single bad row in a vector stops the whole call, and NA counts as bad.
  expect_error(bipois_lpmf_r(3L, 4L, 5, 3, 4, y1_obs = c(1L, 0L, 7L)),
               "found 7", fixed = TRUE)
  expect_error(bipois_lpmf_r(3L, 4L, 5, 3, 4, y1_obs = NA),
               "found NA", fixed = TRUE)
  expect_error(binegbin_lpmf_r(3L, 4L, 5, 3, 4, 2, 1.5, y1_obs = -1L),
               "found -1", fixed = TRUE)

  # The valid values still evaluate, as integers or as doubles.
  expect_true(all(is.finite(bipois_lpmf_r(3L, 4L, 5, 3, 4, y1_obs = c(0, 1)))))
  expect_true(all(is.finite(
    binegbin_lpmf_r(3L, 4L, 5, 3, 4, 2, 1.5, y1_obs = c(0L, 1L))
  )))
})

# ---------------------------------------------------------------------------
# Stan tests -- opt-in via NOT_CRAN=true
# ---------------------------------------------------------------------------

stan_code <- paste0(
  "functions {\n", bipois_stan_funs, binegbin_stan_funs, "}\nmodel {}\n"
)

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

test_that("the exposed Stan functions reject an indicator of 2 and accept 0 and 1", {
  skip_if_not(stan_ready, "rstan unavailable or Stan compilation failed")

  msg <- "y1_obs must be 0 or 1; found 2"
  expect_error(bipois_lpmf(3L, 5, 3, 4, 4L, 2L), msg, fixed = TRUE)
  expect_error(binegbin_lpmf(3L, 5, 3, 4, 2, 1.5, 1.5, 4L, 2L),
               msg, fixed = TRUE)
  expect_error(bipois_lpmf(3L, 5, 3, 4, 4L, -1L), "found -1", fixed = TRUE)

  for (flag in c(0L, 1L)) {
    expect_true(is.finite(bipois_lpmf(3L, 5, 3, 4, 4L, flag)))
    expect_true(is.finite(binegbin_lpmf(3L, 5, 3, 4, 2, 1.5, 1.5, 4L, flag)))
  }
})
