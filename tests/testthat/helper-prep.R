# Builds a minimal synthetic brms "prep"-like list sufficient for
# posterior_predict_<family>()/posterior_epred_<family>() to run without a
# full brmsfit -- avoids the Stan-compilation dependency the rest of this
# package's Stan-side tests already gate behind stan_ready.
# Confirmed sufficient directly against the installed brms package: with
# class(prep) <- "brmsprep", brms::is.brmsprep() (inherits(x, "brmsprep"))
# and brms::get_dpar() (which reads prep$dpars[[name]] and, for a plain
# numeric matrix, slices column i via slice_col()) both work exactly as
# needed against this hand-built structure -- no lazy-evaluation or S3
# dispatch beyond class(prep) <- "brmsprep" is required.
#
# `dpars` is a named list; each entry may be a full `ndraws x nobs` matrix,
# or a shorter vector/scalar recycled across all observations (a single
# dpar value shared by every observation, the common case in these tests).
# `lb`/`ub` are omitted entirely (not merely set to -Inf/Inf) when NULL, to
# match brms's own behaviour for untruncated formulas (confirmed via
# brms::make_standata(): prep$data$lb/ub are absent, not -Inf/Inf-filled,
# when no trunc()/resp_trunc() is used).
make_synthetic_prep <- function(dpars, Y, lb = NULL, ub = NULL,
                                vint1 = NULL, vint2 = NULL) {
  nobs <- length(Y)
  ndraws <- unique(vapply(dpars, function(x) {
    if (is.matrix(x)) nrow(x) else length(x)
  }, numeric(1)))
  ndraws <- ndraws[ndraws > 1]
  ndraws <- if (length(ndraws) == 0) 1 else unique(ndraws)
  stopifnot(length(ndraws) == 1)

  dpars_mat <- lapply(dpars, function(x) {
    if (is.matrix(x)) return(x)
    matrix(x, nrow = ndraws, ncol = nobs)
  })

  data <- list(Y = Y)
  if (!is.null(lb)) data$lb <- lb
  if (!is.null(ub)) data$ub <- ub
  # Supplementary integer data threaded through brms's vint() addition term;
  # the joint families read these as prep$data$vint1/vint2 (e.g. y2 and the
  # y1_obs branch flag for binegbin_partialobs), row-indexed exactly like Y.
  if (!is.null(vint1)) data$vint1 <- vint1
  if (!is.null(vint2)) data$vint2 <- vint2

  structure(
    list(dpars = dpars_mat, data = data, ndraws = ndraws, nobs = nobs),
    class = "brmsprep"
  )
}
