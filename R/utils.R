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
  .get_dpar_any(prep, c(new, old), i = i)
}

# --------------------------------------------------------------------------
# Reading a dpar under any of several accepted spellings
# --------------------------------------------------------------------------
#
# The generalisation of .get_rate() to an ordered list of candidate names,
# needed once 0.8.0 split binegbin_joint's single excess dispersion `shapex`
# into the per-margin pair `shapexone`/`shapextwo` (see NEWS). A stored fit
# may spell that dispersion any of three ways, and which one it uses is a
# property of the fit, not of the attached package:
#
#   shapexone <- "shapexone" (0.8.0+) | "shapexem" (project-local ax family)
#                | "shapex"   (<= 0.7.0 symmetric fits, where one dispersion
#                              served both margins)
#
# Falling back to `shapex` LAST is what makes the symmetric special case work
# without a shim: a five-dpar fit resolves both shapexone and shapextwo to its
# single `shapex`, which is precisely the constraint shapexone == shapextwo
# that the five-dpar family imposed.
#
# The rationale for resolving names against the fit rather than rewriting prep
# or aliasing at the family level is given at .get_rate() above and applies
# unchanged here: brms resolves a dpar's LINK by looking its name up in
# prep$family, so the name handed to brms::get_dpar() must be one the fit's own
# family object declares.
.get_dpar_any <- function(prep, candidates, i = NULL) {
  for (nm in candidates) {
    if (!is.null(prep$dpars[[nm]])) {
      return(brms::get_dpar(prep, nm, i = i))
    }
  }
  stop(
    "This fit has none of the accepted spellings of this distributional ",
    "parameter (", paste0("`", candidates, "`", collapse = ", "),
    "). Available dpars: ", paste(names(prep$dpars), collapse = ", "), ".",
    call. = FALSE
  )
}
