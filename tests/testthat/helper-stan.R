# ==========================================================================
# Stan toolchain check
# ==========================================================================
#
# An installed rstan does not imply a working C++ toolchain, so
# skip_if_not_installed("brms") does not guard a test that fits a model. The
# check below compiles an empty Stan program once per test run and caches the
# result.

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
# Opt-in switch for Stan compilation
# --------------------------------------------------------------------------
#
# Several test files compile a model at file level, outside any test_that(),
# so that rstan::expose_stan_functions() can expose an lpmf for comparison
# with the R reference. No skip_*() reaches such a block, so the block is
# guarded by stan_tests_enabled(), which reads the same NOT_CRAN switch as
# skip_on_cran(). With NOT_CRAN unset, every Stan block is skipped without
# compiling. Rstan declares BH in LinkingTo, so a CI job that installs hard
# dependencies only can still compile Stan.
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
