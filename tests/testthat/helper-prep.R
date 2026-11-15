# Builds a minimal list of class "brmsprep" on which the log_lik_*(),
# posterior_predict_*() and posterior_epred_*() methods run without a fitted
# model, and therefore without Stan compilation. brms::is.brmsprep() and
# brms::get_dpar() require nothing beyond the class and the dpars matrices.
#
# `dpars` is a named list; each entry is an `ndraws x nobs` matrix, or a
# vector or scalar recycled across all observations. `lb` and `ub` are omitted
# when NULL, as in the prep that brms builds for an untruncated formula.
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
  # Supplementary integer data from the brms vint() addition term: vint1 is
  # y2, and vint2 is the observation indicator y1_obs of a _partialobs() fit.
  if (!is.null(vint1)) data$vint1 <- vint1
  if (!is.null(vint2)) data$vint2 <- vint2

  structure(
    list(dpars = dpars_mat, data = data, ndraws = ndraws, nobs = nobs),
    class = "brmsprep"
  )
}
