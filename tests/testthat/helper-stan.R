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
