# ==========================================================================
# Recovery checks: smoke gates vs calibration
# ==========================================================================
#
# WHY THIS EXISTS. The recovery tests used to assert that the true value falls
# inside a single fit's 90% credible interval. For a CORRECT model with a
# well-calibrated posterior that assertion is a Bernoulli(0.9) draw, so it
# fails 10% of the time by construction -- and across the ~18 such assertions
# in this suite roughly 1.8 spurious failures per run were expected. Two were
# observed on every run -- the binegbin recovery assertions on b_lamx_Intercept
# and b_shapex_Intercept -- identically before and after the 0.7.0 dpar rename,
# and went unnoticed from the day they were written because skip_on_cran() keeps
# the fitting tests out of R CMD check and nothing in CI ever set NOT_CRAN.
#
# The two jobs were conflated, so they are now separated:
#
#   recovery_ok()        -- a SMOKE gate. Wide interval (99% by default). Cheap,
#                           runs with the rest of the fitting tests, and answers
#                           "is this fit grossly wrong?". It is deliberately NOT
#                           a calibration statement.
#
#   coverage_recovery()  -- the CALIBRATION instrument. Refits on fresh
#                           simulated data R times and measures how often the
#                           truth lands inside the nominal interval, which is
#                           what a recovery claim actually asserts. Gated behind
#                           an explicit environment variable because it costs
#                           minutes, not seconds.

# --------------------------------------------------------------------------
# Smoke gate
# --------------------------------------------------------------------------

# Shared replacement for the three near-identical local check_recovery()
# definitions that used to live in test-bipois.R, test-binegbin.R and
# test-binegbin_joint.R. `level` is wide on purpose: see above.
recovery_ok <- function(draws, true_val, draws_col, level = 0.99) {
  stopifnot(draws_col %in% names(draws))
  a <- (1 - level) / 2
  q <- stats::quantile(draws[[draws_col]], c(a, 1 - a), names = FALSE)
  true_val >= q[[1]] && true_val <= q[[2]]
}

# --------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------

# The full coverage assessment is minutes of sampling per family. Opt in with
#
#   PAIREDCOUNTBRMS_COVERAGE=true
#
# Run it before a release, or when changing a likelihood, a link, or a
# parameterisation -- not on every save. Using a named variable rather than
# leaning on skip_on_cran() means "this did not run" is a deliberate,
# inspectable state instead of the default nobody notices.
coverage_enabled <- function() {
  isTRUE(as.logical(Sys.getenv("PAIREDCOUNTBRMS_COVERAGE", "false")))
}

skip_unless_coverage <- function() {
  testthat::skip_if_not(
    coverage_enabled(),
    "set PAIREDCOUNTBRMS_COVERAGE=true to run the coverage assessment"
  )
}

# Smallest number of covering replicates out of R that a correctly-calibrated
# model clears with probability >= 1 - alpha. Derived from the Binomial(R,
# level) null rather than picked by eye, so the false-failure rate is a stated
# quantity: at R = 10, level = 0.9, alpha = 0.01 this returns 6, and
# P(X < 6) = P(X <= 5) = 0.0016 -- a 0.16% false-failure rate, against 10% per
# assertion for the single-fit check it replaces.
#
# The floor stays sensitive to real breakage: a mis-specified likelihood or a
# swapped parameter drives coverage towards 0, nowhere near 6/10.
coverage_floor <- function(R, level = 0.9, alpha = 0.01) {
  stats::qbinom(alpha, size = R, prob = level)
}

# Refit `fit` on R fresh datasets from `sim` and return, per parameter, the
# proportion of replicates whose `level` interval contained the truth.
#
# `fit`    a brmsfit to reuse -- its COMPILED Stan model is what makes this
#          affordable. update(recompile = FALSE) re-runs sampling only;
#          measured on this package's binegbin model, compiling costs ~66 s and
#          each refit ~15 s, so R = 10 is ~3.5 min rather than ~11. stanvars
#          survive the update (verified), which is why a custom family can be
#          driven this way at all.
# `sim`    function(i) -> data.frame with the same columns as the original data.
# `truths` named numeric; names are draws-column names, values the true values
#          ON THE SCALE OF THOSE COLUMNS (i.e. already log-transformed for
#          log-linked dpars).
# Extra arguments in `...` are passed to update(), so the replicates can run
# shorter chains than the original fit: coverage needs the 5%/95% quantiles of
# each marginal, not high-precision means, so a couple of thousand draws is
# ample and costs a fraction of the parent fit.
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
