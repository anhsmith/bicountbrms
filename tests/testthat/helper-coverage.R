# ==========================================================================
# Recovery checks: smoke gate and calibration
# ==========================================================================
#
# A single-fit check that a 90% interval contains the true value fails with
# probability 0.1 for a correct model. recovery_ok() is therefore a smoke gate
# with a 99% interval, which detects a grossly wrong fit. coverage_recovery()
# is the calibration check: it refits on R fresh simulated datasets and
# reports the proportion of replicates in which the interval contains the
# true value.
#
# Each fitting test asserts max Rhat < 1.02 rather than the 1.01 recommended
# for inference (Vehtari et al. 2021). On the bipois_partialobs() fit, with
# the data held fixed, max Rhat ranged from 1.0023 to 1.0079 across Stan seeds
# and reached 1.012 on Linux CI, with zero divergent transitions in every run.
# The gate of 1.02 lies below the value of 1.05 at which brms warns, and
# divergent transitions are asserted separately at exactly zero.

# --------------------------------------------------------------------------
# Smoke gate
# --------------------------------------------------------------------------

# TRUE when the central `level` interval of `draws_col` contains `true_val`.
recovery_ok <- function(draws, true_val, draws_col, level = 0.99) {
  stopifnot(draws_col %in% names(draws))
  a <- (1 - level) / 2
  q <- stats::quantile(draws[[draws_col]], c(a, 1 - a), names = FALSE)
  true_val >= q[[1]] && true_val <= q[[2]]
}

# --------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------

# The coverage assessment takes several minutes of sampling per family and
# runs only when BICOUNTBRMS_COVERAGE=true.
coverage_enabled <- function() {
  isTRUE(as.logical(Sys.getenv("BICOUNTBRMS_COVERAGE", "false")))
}

skip_unless_coverage <- function() {
  testthat::skip_if_not(
    coverage_enabled(),
    "set BICOUNTBRMS_COVERAGE=true to run the coverage assessment"
  )
}

# Smallest number of covering replicates out of R that a correctly-calibrated
# model clears with probability >= 1 - alpha, from the Binomial(R, level)
# distribution. At R = 10, level = 0.9 and alpha = 0.01, the floor is 6, and
# P(X <= 5) = 0.0016. A mis-specified likelihood or a swapped parameter drives
# coverage towards 0.
coverage_floor <- function(R, level = 0.9, alpha = 0.01) {
  stats::qbinom(alpha, size = R, prob = level)
}

# Refit `fit` on R fresh datasets from `sim` and return, per parameter, the
# proportion of replicates whose `level` interval contained the truth.
#
# `fit`    a brmsfit to reuse. update(recompile = FALSE) re-runs sampling
#          without recompiling, and the stanvars of the custom family are kept.
# `sim`    function(i) -> data.frame with the same columns as the original data.
# `truths` named numeric; names are draws-column names, and values are the
#          true values on the scale of those columns (log scale for log-linked
#          dpars).
# Extra arguments in `...` are passed to update(), so the replicates can run
# shorter chains than the original fit; coverage requires only the 5% and 95%
# quantiles of each marginal posterior.
coverage_recovery <- function(fit, sim, truths, R = 10L, level = 0.9,
                              verbose = TRUE, ...) {
  a <- (1 - level) / 2
  hits <- matrix(
    NA, nrow = R, ncol = length(truths),
    dimnames = list(NULL, names(truths))
  )

  for (r in seq_len(R)) {
    f <- suppressMessages(suppressWarnings(
      stats::update(fit, newdata = sim(r), refresh = 0, recompile = FALSE, ...)
    ))
    d <- as.data.frame(f)
    for (p in names(truths)) {
      q <- stats::quantile(d[[p]], c(a, 1 - a), names = FALSE)
      hits[r, p] <- truths[[p]] >= q[[1]] && truths[[p]] <= q[[2]]
    }
    if (verbose) {
      message(sprintf("  coverage replicate %d/%d done", r, R))
    }
  }

  colMeans(hits)
}
