# ==========================================================================
# Deprecated names retained from 0.4.0-0.8.0
#
# 0.9.0 renamed the censoring-aware family binegbin_joint -> binegbin_cens.
# The old suffix said the family models the pair jointly, which is true of all
# four families in this package and so distinguished nothing; `_cens` names the
# property that actually selects it, that one margin may be unobserved.
#
# TWO SEPARATE COMPATIBILITY PROBLEMS, WITH DIFFERENT SOLUTIONS.
#
# Source compatibility -- code that calls binegbin_joint() -- is handled by the
# constructor and stanvars shims below. They return the NEW family, so a fit
# made through the deprecated name is in every respect a binegbin_cens fit and
# carries no legacy naming forward.
#
# Fit compatibility is the harder one and is why the three post-processing
# shims exist. A brmsfit stores its own family object, so a fit made before
# 0.9.0 has family$name == "binegbin_joint" permanently. brms resolves
# post-processing methods by name off the search path at call time --
# log_lik_<name>, posterior_predict_<name>, posterior_epred_<name> -- built
# from the STORED name, not from whatever the attached package now declares.
# Without these three, loo(), posterior_predict() and posterior_epred() on a
# stored fit fail with "could not find function", and no argument to those
# calls can route around it. Every dpar name is unchanged by the rename, so
# forwarding is a straight hand-off.
#
# The post-processing shims deliberately do NOT warn. brms calls log_lik_ and
# posterior_predict_ once per observation, so a warning there would fire
# thousands of times per loo() call, and it would be unactionable in any case:
# the name is baked into a fit the user cannot edit. The constructor and
# stanvars shims warn, because those are the calls a user can actually change.
#
# Removal is deferred to the next major version. bipois_cens() gets no shims:
# it was added in this same release under its final name, so no fit and no
# source can refer to it as bipois_joint().

#' Deprecated names for the censoring-aware Negative-Binomial family
#'
#' @description
#' `binegbin_joint()` and `binegbin_joint_stanvars()` were renamed to
#' [binegbin_cens()] and [binegbin_cens_stanvars()] in 0.9.0. The old names
#' still work and return the new objects, with a deprecation warning.
#'
#' `log_lik_binegbin_joint()`, `posterior_predict_binegbin_joint()` and
#' `posterior_epred_binegbin_joint()` exist so that fits made before 0.9.0
#' keep working. Such a fit stores `family$name == "binegbin_joint"`, and brms
#' builds its post-processing method names from that stored name, so these
#' three are reached by `loo()`, `posterior_predict()` and `posterior_epred()`
#' without the user naming them. They forward silently -- the name is a
#' property of the stored fit, not of anything the caller can change, and brms
#' calls the first two once per observation. **No refitting is required.**
#'
#' New code should use [binegbin_cens()]. These names will be removed in the
#' next major version.
#'
#' @param i,prep Passed through unchanged to the corresponding
#'   `binegbin_cens` method.
#' @param ... Passed through unchanged.
#'
#' @return As the corresponding [binegbin_cens()] function.
#' @name binegbin_joint-deprecated
#' @keywords internal
NULL

#' @rdname binegbin_joint-deprecated
#' @export
binegbin_joint <- function() {
  .Deprecated("binegbin_cens", package = "bicountbrms")
  binegbin_cens()
}

#' @rdname binegbin_joint-deprecated
#' @export
binegbin_joint_stanvars <- function() {
  .Deprecated("binegbin_cens_stanvars", package = "bicountbrms")
  binegbin_cens_stanvars()
}

#' @rdname binegbin_joint-deprecated
#' @export
log_lik_binegbin_joint <- function(i, prep) {
  log_lik_binegbin_cens(i, prep)
}

#' @rdname binegbin_joint-deprecated
#' @export
posterior_predict_binegbin_joint <- function(i, prep, ...) {
  posterior_predict_binegbin_cens(i, prep, ...)
}

#' @rdname binegbin_joint-deprecated
#' @export
posterior_epred_binegbin_joint <- function(prep) {
  posterior_epred_binegbin_cens(prep)
}
