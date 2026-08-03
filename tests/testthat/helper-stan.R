# ==========================================================================
# Is there a working Stan toolchain?
# ==========================================================================
#
# `skip_if_not_installed("brms")` is NOT sufficient to guard a test that fits a
# model. brms imports rstan, so both are present whenever brms is, but rstan
# still cannot COMPILE a model without the Boost headers (package BH) and a C++
# toolchain. A CI runner installing hard dependencies only gets brms and rstan
# and no BH, so `brm()` reaches `rstan::stan_model()` and fails with
# "Boost not found" rather than skipping.
#
# That is exactly what happened on the first CI run of this package: the
# difference-family fitting tests skipped cleanly because they guard on a
# file-local `stan_ready`, while the joint-family ones and test-recovery.R had
# only skip_on_cran() + skip_if_not_installed("brms") and errored. Note also
# that r-lib/actions/check-r-package sets NOT_CRAN=true, so skip_on_cran() does
# NOT fire in that workflow.
#
# The check below compiles the smallest possible Stan program, which is the only
# way to establish that compilation actually works rather than that rstan is
# merely installed. The result is cached, so the cost is paid once per test run
# rather than once per file.

.stan_toolchain_ok <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ok <- FALSE
    if (requireNamespace("rstan", quietly = TRUE)) {
      ok <- tryCatch(
        {
          suppressMessages(suppressWarnings(
            rstan::stan_model(model_code = "model {}")
          ))
          TRUE
        },
        error = function(e) FALSE
      )
    }
    cached <<- ok
    ok
  }
})

# Guard for any test that calls brm(). Use alongside skip_on_cran(); neither
# implies the other.
skip_if_no_stan <- function() {
  testthat::skip_if_not(
    .stan_toolchain_ok(),
    "no working Stan toolchain (rstan cannot compile a model)"
  )
}

# --------------------------------------------------------------------------
# Whether to ATTEMPT Stan compilation at all
# --------------------------------------------------------------------------
#
# The two guards above decide whether a Stan test can run. This one decides
# whether it is worth asking, and it exists because the expensive part of this
# suite is not the sampling -- it is compilation.
#
# Several test files compile a model at FILE level, outside any test_that(), so
# that rstan::expose_stan_functions() can expose an lpmf for comparison against
# the R reference. Being outside a test, no skip_*() can reach them: they run
# whenever the file is sourced, and they cost minutes each.
#
# That was tolerable while CI could not compile at all. It cannot be relied on:
# `dependencies: '"hard"'` installs brms -> rstan, and rstan declares
# `LinkingTo: BH, StanHeaders, ...` -- LinkingTo IS a hard dependency, so the
# Boost headers arrive with it and rstan compiles fine. The fast CI job was
# therefore compiling every model and running every fit, taking ~40 minutes
# against the ~5 its comments claimed.
#
# So compilation is now opt-in, on the same switch that governs the fits
# (NOT_CRAN), which keeps "did the Stan tests run?" a single question with a
# single answer rather than two that can disagree. Set NOT_CRAN=true to compile;
# leave it unset and every Stan block skips without attempting a compile, so a
# push gets the analytic and R-side signal in a couple of minutes.
#
# Use as the condition on a file-level compile block:
#
#   stan_ready <- FALSE
#   if (stan_tests_enabled() && requireNamespace("rstan", quietly = TRUE)) {
#     tryCatch({ ...compile and expose... ; stan_ready <- TRUE },
#              error = function(e) NULL)
#   }
stan_tests_enabled <- function() {
  isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))
}
