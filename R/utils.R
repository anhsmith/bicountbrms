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
# Every rate read in every family goes through here -- bipois.R,
# bipois_cens.R, binegbin.R and binegbin_cens.R --
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
# needed once 0.8.0 split binegbin_cens's single excess dispersion `shapex`
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

# --------------------------------------------------------------------------
# E[N_shared | y2] for the negative-binomial joint families
# --------------------------------------------------------------------------
#
# Every joint family's conditional expectation has the same shape,
#
#   E[y1 | y2] = E[N_shared | y2] + lambdaone,
#
# and the families differ only in the first term. For bipois()/bipois_cens()
# it is closed form -- conditioning a sum of independent Poissons on its total
# gives a Binomial, so E[N_shared | y2] = y2 * mu/(mu + lambdatwo) -- and those
# families compute it inline. A sum of independent negative binomials admits no
# such shortcut, so binegbin()/binegbin_cens() need the conditional evaluated
# directly:
#
#   P(N_shared = k | y2) proportional to
#     NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex2),  k = 0..y2
#
#   E[N_shared | y2] = sum_k k P(N_shared = k | y2)
#
# These are the same weights posterior_predict_binegbin() and
# posterior_predict_binegbin_cens() already build; this helper sums them
# instead of sampling from them. That makes the epred exact rather than an
# approximation, and makes it consistent with the predictions by construction
# -- the two now differ only by Monte Carlo error, which is what a test can
# check.
#
# Until 0.9.0 posterior_epred_binegbin() substituted the MARGINAL shared
# fraction mu/(mu + lambdatwo) -- the bipois answer -- for the conditional one,
# and posterior_epred_binegbin_cens() did not exist. The substitution is exact
# only in the Poisson limit and biased otherwise, in the direction set by which
# component carries more dispersion. See NEWS for the size of the change.
#
# `shapex2` is the dispersion of the SECOND margin's excess component
# (`shapex` for binegbin(), `shapextwo` for binegbin_cens()) -- the y2
# marginal is what is being conditioned on, so the first margin's dispersion
# does not enter.
#
# Arguments are ndraws x nobs matrices except `y2`, a length-nobs integer
# vector. Returns an ndraws x nobs matrix. The loop is over observations and
# over the support 0..y2, vectorised across draws within each, so its cost
# matches posterior_predict's; both are paid post-hoc, never inside the
# sampler.
.e_shared_given_y2_nb <- function(mu, lambdatwo, shapes, shapex2, y2) {
  out <- matrix(0, nrow = nrow(mu), ncol = ncol(mu))
  for (j in seq_along(y2)) {
    yj <- y2[[j]]
    if (yj == 0L) next            # N_shared | y2 = 0 is degenerate at 0
    k <- 0:yj
    # log-weights: ndraws x (yj + 1), one column per k.
    lw <- vapply(k, function(kk) {
      stats::dnbinom(kk,      size = shapes[, j],  mu = mu[, j],        log = TRUE) +
        stats::dnbinom(yj - kk, size = shapex2[, j], mu = lambdatwo[, j], log = TRUE)
    }, numeric(nrow(mu)))
    if (!is.matrix(lw)) lw <- matrix(lw, nrow = nrow(mu))
    w <- exp(lw - apply(lw, 1, max))
    out[, j] <- (w %*% k) / rowSums(w)
  }
  out
}
