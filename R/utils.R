# ==========================================================================
# Internal helpers
# ==========================================================================

# --------------------------------------------------------------------------
# Reading a rate dpar under its current or pre-0.7.0 name
# --------------------------------------------------------------------------
#
# Fits made before 0.7.0 declare the exclusive rates as lambdaem and lambdalb
# rather than lambdaone and lambdatwo. A brmsfit stores its own family object,
# and prepare_predictions() builds prep$dpars from that stored family, so the
# prep of such a fit has prep$dpars$lambdaem and no prep$dpars$lambdaone.
#
# The name passed to brms::get_dpar() must be declared by the stored family,
# not merely present in prep$dpars, because brms:::apply_dpar_inv_link()
# finds the link of a dpar by looking its name up in prep$family. For a name
# that is absent, brms::get_dpar() stops with the bare message
# `!is.null(x) is not TRUE`. .get_rate() therefore passes brms the first name
# that the stored family declares. Every post-processing method in bipois.R,
# binegbin.R, and predict.R reads the exclusive rates through .get_rate().
.get_rate <- function(prep, new, old, i = NULL) {
  .get_dpar_any(prep, c(new, old), i = i)
}

# --------------------------------------------------------------------------
# Reading a dpar under any of several accepted spellings
# --------------------------------------------------------------------------
#
# .get_rate() generalised to an ordered list of candidate names, for the
# dispersions of the exclusive components. A stored binegbin fit declares
# the dispersion of the first exclusive component under one of three names:
#
#   "shapexone"  current fits
#   "shapexem"   fits from an earlier asymmetric family that declared
#                shapexem and shapexlb
#   "shapex"     binegbin fits made before 0.10.0, with one dispersion for
#                both exclusive components
#
# `shapex` is the last candidate in both .SHAPEXONE_NAMES and
# .SHAPEXTWO_NAMES (binegbin.R), so a five-dpar fit reads its single `shapex`
# as both shapexone and shapextwo. That is the constraint
# shapexone == shapextwo under which the five-dpar fit was made. The reason
# for passing brms a name declared by the stored family is given at
# .get_rate().
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
# For both families, the conditional expectation of the first count is
#
#   E[y1 | y2] = E[N_shared | y2] + lambdaone,
#
# and the two families differ only in E[N_shared | y2]. For bipois, conditioning
# a sum of independent Poissons on its total gives a Binomial, so
# E[N_shared | y2] = y2 * mu/(mu + lambdatwo), computed inline in predict.R.
# A sum of independent negative binomials has no such closed form, so for
# binegbin the conditional distribution is evaluated directly:
#
#   P(N_shared = k | y2) proportional to
#     NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex2),  k = 0..y2
#
#   E[N_shared | y2] = sum_k k P(N_shared = k | y2)
#
# The weights P(N_shared = k | y2) are the same weights that
# .impute_draw_binegbin() in predict.R samples from; this helper sums them
# instead. Summing the weights makes posterior_epred_y1_given_y2() exact
# rather than an approximation, and consistent with
# posterior_predict_y1_given_y2() by construction. test-impute.R tests that the
# expectation equals the mean of the draws to within Monte Carlo error.
#
# The marginal shared fraction mu/(mu + lambdatwo) gives the right answer only
# in the Poisson limit, and must not replace the conditional weights above.
#
# `shapex2` is the dispersion of the second exclusive component (`shapextwo`,
# or the single `shapex` of a pre-0.10.0 fit). The expectation is conditional
# on y2, so the dispersion of the first exclusive component does not enter.
#
# Arguments are ndraws x nobs matrices except `y2`, a length-nobs integer
# vector. Returns an ndraws x nobs matrix. The loops run over observations and
# over the support 0..y2, and are vectorised across draws.
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

# --------------------------------------------------------------------------
# The observation flag, for a prep built by either constructor
# --------------------------------------------------------------------------

# Both constructors of a family return the same family name. bipois() and
# binegbin() declare vars = c("vint1[n]", "1"), so the flag reaches Stan as a
# literal and never enters the data block. A prep built from such a fit has no
# vint2, and every row is paired. bipois_partialobs() and
# binegbin_partialobs() declare vars = c("vint1[n]", "vint2[n]"), and vint2
# holds the flag.
#
# Only log_lik_* calls this. posterior_predict_* and posterior_epred_* return
# the marginal distribution and mean of y1, which do not depend on the flag.
# The _y1_given_y2 functions in predict.R read the flag through .paired_rows(),
# to return the recorded y1 on paired rows when keep_observed = TRUE.
.y1_obs_at <- function(prep, i) {
  if (is.null(prep$data$vint2)) 1L else prep$data$vint2[i]
}

# The R references stop on a flag other than 0 or 1, as the Stan lpmf rejects
# one. brms binds vint() terms by position, so vint(y1_obs, y2) passes the flag
# as the second count and the second count as the flag. log_lik() with newdata
# reaches the R references without passing through Stan, and without this
# check both references would evaluate any flag other than 1 on the unpaired
# branch and return a finite log-density. NA fails the check too, because %in%
# returns FALSE for it.
.check_y1_obs <- function(y1_obs) {
  bad <- !(y1_obs %in% c(0, 1))
  if (any(bad)) {
    stop(
      "y1_obs must be 0 or 1; found ", y1_obs[bad][1],
      ". Check the order: vint(y2, y1_obs).",
      call. = FALSE
    )
  }
  invisible(y1_obs)
}
