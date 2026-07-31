# ==========================================================================
# (M, f, delta) <-> native dpar coordinates for the joint bivariate families
#
# binegbin()/bipois() are parameterised by three rates -- mu (shared),
# lambdaone and lambdatwo (the two source-specific excesses) -- because that
# is what the trivariate-reduction likelihood consumes directly. Those three
# are correlated in use: raising the overall catch level moves all three at
# once, so none of them is individually interpretable as "how much was
# caught", "how much did the two sources agree", or "which source ran high".
#
# The (M, f, delta) coordinates separate exactly those three questions:
#
#   M     = mu + (lambdaone + lambdatwo)/2  overall level (midpoint of the two
#                                           sources' expectations)
#   f     = mu / M                          congruence: the share of M that
#                                           both sources saw
#   delta = 0.5 * log(lambdaone/lambdatwo)  source bias, on a log-ratio scale
#
# with inverse
#
#   mu        = M f
#   lambdaone = M (1 - f) (1 + tanh delta)
#   lambdatwo = M (1 - f) (1 - tanh delta)
#
# The bounded bias beta = tanh(delta) in [-1, 1] is often the more convenient
# dial: beta = (lambdaone - lambdatwo)/(lambdaone + lambdatwo) reads directly
# as the fractional imbalance of excess between the two sources, and
# beta = +/-1 is the limit where one source never records an unshared event.
#
# The map is a bijection on the interior. Because the average of the two excess
# rates is M(1-f) for ANY delta, M stays pinned to the midpoint whatever the
# bias -- M is a midpoint of what the two sources REPORT, not a property of the
# underlying process.
#
# DISPERSION. The dpars shapes/shapex are NB2 phi (Stan neg_binomial_2, R
# dnbinom size): variance = m + m^2/phi, so LARGER phi means LESS
# overdispersion and Poisson is the phi -> Inf limit. The SD-scale kappa used
# alongside (M, f, delta) inverts that into a dial that increases with
# overdispersion and reaches Poisson at a finite zero:
#
#   shapes = 1/kappas^2      kappas = 1/sqrt(shapes)
#
# Note the direction reversal: raising kappa LOWERS shape.
#
# These are pure coordinate transforms -- they fit nothing. To fit in (M, f,
# delta) coordinates, supply them through a non-linear formula; see the
# examples on binegbin_mfd_to_dpars(). Use these helpers to set up simulations
# from interpretable values, and to read fitted dpars back into interpretable
# ones.

# --------------------------------------------------------------------------

#' Convert (M, f, delta) coordinates to native binegbin/bipois dpars
#'
#' @description
#' Maps the interpretable coordinates -- overall level `M`, congruence `f`, and
#' source bias `delta` -- onto the rate dpars that [binegbin()] and [bipois()]
#' actually take (`mu`, `lambdaone`, `lambdatwo`), optionally converting
#' SD-scale dispersions `kappas`/`kappax` to the `shapes`/`shapex` dpars.
#'
#' [binegbin_dpars_to_mfd()] is the exact inverse.
#'
#' @param M Overall level: `mu + (lambdaone + lambdatwo)/2`. Non-negative.
#' @param f Congruence, the share of `M` that both sources saw: `mu / M`. In
#'   `[0, 1]`. `f = 1` means perfect agreement (both excesses vanish); `f = 0`
#'   means no shared component at all.
#' @param delta Source bias on the log-ratio scale,
#'   `0.5 * log(lambdaone/lambdatwo)`. `0` is unbiased. `+/-Inf` is permitted
#'   and gives the limit where one excess rate is zero.
#' @param kappas,kappax Optional SD-scale dispersions for the shared and excess
#'   components. `0` is the Poisson limit. If supplied, the returned list gains
#'   `shapes`/`shapex` (`= 1/kappa^2`, so `kappa = 0` gives `Inf`).
#'
#' @details
#' Arguments are recycled to a common length, so this vectorises over posterior
#' draws.
#'
#' **Boundary behaviour.** At `f = 1` both excess rates are exactly `0`
#' regardless of `delta` -- the bias becomes unidentifiable, which
#' [binegbin_dpars_to_mfd()] reports back as `NA`. This direction is always
#' well defined; only the inverse degenerates.
#'
#' @return A named list of `mu`, `lambdaone`, `lambdatwo`, plus `shapes` and
#'   `shapex` when `kappas`/`kappax` are supplied.
#'
#' @examples
#' # A moderately congruent pair, source 1 running high
#' binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#'
#' # Perfect congruence: both excesses vanish
#' binegbin_mfd_to_dpars(M = 12, f = 1, delta = 0.5)
#'
#' # To FIT in these coordinates, pass them through a non-linear formula
#' # (all five dpars are log-linked, so the link supplies the exp()):
#' #   bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#' #     nlf(mu        ~ eta + log_inv_logit(con)) +
#' #     nlf(lambdaone ~ log(2) + eta + log_inv_logit(-con) +
#' #                     log_inv_logit(2 * methd)) +
#' #     nlf(lambdatwo ~ log(2) + eta + log_inv_logit(-con) +
#' #                     log_inv_logit(-2 * methd)) +
#' #     lf(eta ~ 1, con ~ 1, methd ~ 1)
#' # where eta = log M, con = logit f, methd = delta.
#'
#' @seealso [binegbin_dpars_to_mfd()], [binegbin()], [bipois()]
#' @export
binegbin_mfd_to_dpars <- function(M, f, delta = 0, kappas = NULL, kappax = NULL) {
  n <- max(length(M), length(f), length(delta))
  M     <- rep_len(M, n)
  f     <- rep_len(f, n)
  delta <- rep_len(delta, n)

  if (any(M < 0, na.rm = TRUE)) stop("`M` must be non-negative.", call. = FALSE)
  if (any(f < 0 | f > 1, na.rm = TRUE)) stop("`f` must lie in [0, 1].", call. = FALSE)

  excess_mid <- M * (1 - f)          # = (lambdaone + lambdatwo)/2, for any delta
  b <- tanh(delta)                   # bounded bias in [-1, 1]; tanh(+-Inf) = +-1

  out <- list(
    mu        = M * f,
    lambdaone = excess_mid * (1 + b),
    lambdatwo = excess_mid * (1 - b)
  )

  if (!is.null(kappas)) out$shapes <- .kappa_to_shape(rep_len(kappas, n))
  if (!is.null(kappax)) out$shapex <- .kappa_to_shape(rep_len(kappax, n))
  out
}

#' Convert native binegbin/bipois dpars to (M, f, delta) coordinates
#'
#' @description
#' Exact inverse of [binegbin_mfd_to_dpars()]. Reads the rate dpars `mu`,
#' `lambdaone`, `lambdatwo` back into the interpretable overall level `M`,
#' congruence `f`, and source bias `delta`, optionally converting
#' `shapes`/`shapex` back to SD-scale `kappas`/`kappax`.
#'
#' @param mu Shared-component rate.
#' @param lambdaone,lambdatwo The two excess rates.
#' @param shapes,shapex Optional NB2 dispersions. If supplied, the returned list
#'   gains `kappas`/`kappax` (`= 1/sqrt(shape)`, so `shape = Inf` gives `0`).
#'
#' @details
#' Arguments are recycled to a common length, so this vectorises over posterior
#' draws -- e.g. to convert a whole posterior into interpretable coordinates.
#'
#' **Boundary behaviour**, which the forward direction does not have:
#'
#' * `lambdaone == lambdatwo == 0` (perfect congruence, `f = 1`): the bias is
#'   genuinely unidentified -- there is no excess to be biased -- and `delta` is
#'   returned as `NA`, not `0`. Zero would assert an unbiased source, which the
#'   data at that point cannot support.
#' * `M == 0` (nothing anywhere): `f` is undefined and returned as `NA`.
#' * Exactly one excess rate `0`: `delta` is `+/-Inf`, the well-defined limit
#'   where one source never records an unshared event.
#'
#' Because of the first case, round-tripping is exact everywhere except at
#' `f = 1`, where `delta` cannot be recovered. Hold one coordinate system as the
#' source of truth rather than repeatedly converting back and forth.
#'
#' @return A named list of `M`, `f`, `delta`, plus `beta` (the bounded bias
#'   `tanh(delta)`), and `kappas`/`kappax` when `shapes`/`shapex` are supplied.
#'
#' @examples
#' binegbin_dpars_to_mfd(mu = 8.04, lambdaone = 4.75, lambdatwo = 3.17)
#'
#' # Round trip
#' d <- binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#' binegbin_dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)[c("M", "f", "delta")]
#'
#' # Perfect congruence: bias is unidentified, reported as NA
#' binegbin_dpars_to_mfd(mu = 12, lambdaone = 0, lambdatwo = 0)$delta
#'
#' @seealso [binegbin_mfd_to_dpars()], [binegbin()], [bipois()]
#' @export
binegbin_dpars_to_mfd <- function(mu, lambdaone, lambdatwo,
                                  shapes = NULL, shapex = NULL) {
  n <- max(length(mu), length(lambdaone), length(lambdatwo))
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  if (any(c(mu, lambdaone, lambdatwo) < 0, na.rm = TRUE)) {
    stop("Rates must be non-negative.", call. = FALSE)
  }

  excess_sum <- lambdaone + lambdatwo
  M <- mu + excess_sum / 2

  # f undefined when there is nothing at all.
  f <- ifelse(M > 0, mu / M, NA_real_)

  # delta undefined when there is no excess to be biased (f == 1). One rate
  # zero and the other positive is the legitimate +/-Inf limit, which log()
  # produces directly.
  delta <- ifelse(excess_sum > 0, 0.5 * log(lambdaone / lambdatwo), NA_real_)
  beta  <- ifelse(excess_sum > 0, (lambdaone - lambdatwo) / excess_sum, NA_real_)

  out <- list(M = M, f = f, delta = delta, beta = beta)

  if (!is.null(shapes)) out$kappas <- .shape_to_kappa(rep_len(shapes, n))
  if (!is.null(shapex)) out$kappax <- .shape_to_kappa(rep_len(shapex, n))
  out
}

# Internal dispersion conversions. Note the direction reversal: kappa increases
# with overdispersion, shape decreases. kappa = 0 <-> shape = Inf is the
# Poisson limit and is handled exactly rather than by division blowing up.
.kappa_to_shape <- function(kappa) {
  if (any(kappa < 0, na.rm = TRUE)) {
    stop("`kappa` must be non-negative.", call. = FALSE)
  }
  ifelse(kappa == 0, Inf, 1 / kappa^2)
}

.shape_to_kappa <- function(shape) {
  if (any(shape < 0, na.rm = TRUE)) {
    stop("`shape` must be non-negative.", call. = FALSE)
  }
  ifelse(is.infinite(shape), 0, 1 / sqrt(shape))
}
