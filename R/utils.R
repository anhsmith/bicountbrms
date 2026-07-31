# ==========================================================================
# Internal helpers
# ==========================================================================

# --------------------------------------------------------------------------
# Reading a rate dpar across the 0.7.0 rename
# --------------------------------------------------------------------------
#
# 0.7.0 renamed the joint families' two excess-rate dpars lambdaem/lambdalb
# to lambdaone/lambdatwo (see NEWS). Fits made before that carry the OLD
# names, and they carry them in a way post-processing cannot route around
# generically: a brmsfit stores its own family object, so
# prepare_predictions() builds prep$dpars from THAT family, not from
# whatever the currently-attached package declares. An old fit therefore
# arrives here with prep$dpars$lambdaem populated and prep$dpars$lambdaone
# absent.
#
# brms::get_dpar() does `x <- prep$dpars[[dpar]]; stopifnot(!is.null(x))`,
# so asking an old fit for the new name fails on a bare `!is.null(x) is not
# TRUE` with nothing pointing at the cause. Worse, brms:::apply_dpar_inv_link()
# resolves the dpar's LINK by looking the name up in prep$family -- so the
# name passed to get_dpar() must match the fit's own family object, not
# merely be present in prep$dpars. Rewriting prep, or aliasing at the family
# level, would have to satisfy both; reading under the fit's own name
# satisfies both for free.
#
# Hence: resolve the name against the fit, then hand that name to brms.
# Every rate read in bipois.R/binegbin.R/binegbin_joint.R goes through here,
# which is what lets pre-0.7.0 fits keep working with loo(),
# posterior_predict() and log_lik() without refitting.
.get_rate <- function(prep, new, old, i = NULL) {
  nm <- if (!is.null(prep$dpars[[new]])) {
    new
  } else if (!is.null(prep$dpars[[old]])) {
    old
  } else {
    stop(
      "This fit has neither the `", new, "` dpar nor the pre-0.7.0 `", old,
      "`. Available dpars: ", paste(names(prep$dpars), collapse = ", "), ".",
      call. = FALSE
    )
  }
  brms::get_dpar(prep, nm, i = i)
}
