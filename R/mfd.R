# ==========================================================================
# congruence–bias coordinates (M, f, delta) <-> native parameters
#
#   M     = mu + (lambdaone + lambdatwo)/2   midpoint of the two expected counts
#   f     = mu / M                           congruence
#   delta = 0.5 * log(lambdaone/lambdatwo)   source bias, on a log-ratio scale
#
#   mu = M f,   lambdaone = M (1 - f)(1 + tanh delta),
#   lambdatwo = M (1 - f)(1 - tanh delta)
#
# The overdispersion kappa = 1/sqrt(phi) for each NB2 dispersion phi (shapes,
# shapexone, shapextwo); kappa = 0 is the Poisson limit.
#
# The three rates are common to all four constructors. mfd_to_dpars() writes
# only current dpar names; dpars_to_mfd() also reads the single `shapex` of
# binegbin fits made before 0.10.0, as .get_dpar_any() does. Both functions are
# pure transforms. Fitting in these coordinates is shown in
# vignettes/_bicountbrms.Rmd.

# --------------------------------------------------------------------------

#' Convert from the congruence–bias parameterisation to native parameters
#'
#' @description
#' Converts from the congruence–bias parameterisation (the midpoint `M` of the two
#' expected counts, the congruence `f`, and the source bias `delta`) to the
#' native parameters (the rates `mu`, `lambdaone`, and `lambdatwo`), and,
#' optionally, converts overdispersions (`kappas`, `kappax`, `kappaxone`, and
#' `kappaxtwo`) to the native dispersions (`shapes`, `shapexone`, and
#' `shapextwo`). [dpars_to_mfd()] converts in the other direction.
#'
#' @param M Midpoint of the two expected counts:
#'   `mu + (lambdaone + lambdatwo)/2`, equal to `(E[y1] + E[y2])/2`.
#'   Non-negative.
#' @param f Congruence, the mean of the shared component as a proportion of
#'   `M`: `mu / M`. In `[0, 1]`. `f = 1` means perfect congruence (both
#'   exclusive rates are zero); `f = 0` means no shared component.
#' @param delta Source bias on the log-ratio scale,
#'   `0.5 * log(lambdaone/lambdatwo)`. `0` is unbiased. `+/-Inf` is permitted
#'   and gives the limit where one exclusive rate is zero.
#' @param kappas,kappax Optional overdispersions on the scale
#'   `kappa = 1/sqrt(phi)`, where `phi` is the NB2 (inverse) dispersion
#'   parameter; `kappa = 0` is the Poisson limit. `kappas` is the
#'   overdispersion of the shared component. Supplying `kappas` adds
#'   `shapes = 1/kappas^2` to the returned list, so `kappas = 0` gives
#'   `shapes = Inf`. `kappax` is one overdispersion for both exclusive
#'   components. Supplying `kappax` adds `shapexone` and `shapextwo` to the
#'   returned list, both equal to `1/kappax^2`. For [bipois()] and
#'   [bipois_partialobs()], which have no dispersion parameters, leave every
#'   overdispersion argument unset.
#' @param kappaxone,kappaxtwo Optional overdispersions of the two exclusive
#'   components, for the general case in which they differ. Supplying them
#'   adds `shapexone` and `shapextwo` to the returned list. Mutually exclusive
#'   with `kappax`.
#'
#' @details
#' Each argument may be a vector, such as the posterior draws of `M`, `f`, and
#' `delta`. `M`, `f`, and `delta` are recycled to the length of the longest of
#' the three, and any overdispersions to the same length; each element of the
#' returned list is a vector of that length.
#'
#' To fit a model in the congruence–bias parameterisation, use [bf_mfd()], as
#' in `vignette("bicountbrms")`.
#'
#' @return A named list of `mu`, `lambdaone`, and `lambdatwo`, with `shapes`
#'   added when `kappas` is supplied, and `shapexone` and `shapextwo` added
#'   when `kappax`, or `kappaxone` and `kappaxtwo`, are supplied.
#'
#' @examples
#' # A moderately congruent pair, source 1 recording more of the exclusive
#' # events
#' mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#'
#' # Perfect congruence: both exclusive rates are zero
#' mfd_to_dpars(M = 12, f = 1, delta = 0.5)
#'
#' # kappax: one overdispersion for both exclusive components
#' mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4, kappax = 0.9)
#'
#' # kappaxone, kappaxtwo: one overdispersion per exclusive component
#' mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4,
#'              kappaxone = 0.9, kappaxtwo = 0.3)
#'
#' @seealso [dpars_to_mfd()], [bipois()], [bipois_partialobs()],
#'   [binegbin()], [binegbin_partialobs()]
#' @export
mfd_to_dpars <- function(M, f, delta = 0, kappas = NULL, kappax = NULL,
                         kappaxone = NULL, kappaxtwo = NULL) {
  n <- max(length(M), length(f), length(delta))
  M     <- rep_len(M, n)
  f     <- rep_len(f, n)
  delta <- rep_len(delta, n)

  if (any(M < 0, na.rm = TRUE)) stop("`M` must be non-negative.", call. = FALSE)
  if (any(f < 0 | f > 1, na.rm = TRUE)) stop("`f` must lie in [0, 1].", call. = FALSE)
  # `kappax` and the per-margin pair write the SAME two slots -- kappax is the
  # shorthand for giving both margins one value -- so accepting both at once
  # would silently let one overwrite the other.
  if (!is.null(kappax) && (!is.null(kappaxone) || !is.null(kappaxtwo))) {
    stop("Supply either `kappax` (one overdispersion for both exclusive ",
         "components) or `kappaxone`/`kappaxtwo` (a separate overdispersion ",
         "for each exclusive component), not both. ",
         "Both write `shapexone` and `shapextwo`.", call. = FALSE)
  }

  exclusive_mid <- M * (1 - f)  # = (lambdaone + lambdatwo)/2, for any delta
  b <- tanh(delta)              # bounded bias in [-1, 1]; tanh(+-Inf) = +-1

  out <- list(
    mu        = M * f,
    lambdaone = exclusive_mid * (1 + b),
    lambdatwo = exclusive_mid * (1 - b)
  )

  if (!is.null(kappas)) out$shapes <- .kappa_to_shape(rep_len(kappas, n))
  # `kappax` is the tied case, so it writes BOTH per-margin slots. Emitting a
  # `shapex` here would name a dpar no current constructor accepts, and the output
  # of this function is meant to be usable in a brm() call unaltered.
  if (!is.null(kappax)) {
    # Assigned separately, not chained: `a <- b <- value` evaluates
    # right-to-left, so the chained form would create shapextwo first and the
    # returned list would read shapextwo, shapexone, the reverse of the order
    # of the arguments of binegbin_lpmf.
    kx <- .kappa_to_shape(rep_len(kappax, n))
    out$shapexone <- kx
    out$shapextwo <- kx
  }
  if (!is.null(kappaxone)) out$shapexone <- .kappa_to_shape(rep_len(kappaxone, n))
  if (!is.null(kappaxtwo)) out$shapextwo <- .kappa_to_shape(rep_len(kappaxtwo, n))
  out
}

#' Convert from native parameters to the congruence–bias parameterisation
#'
#' @description
#' Converts from the native parameters (the rates `mu`, `lambdaone`, and
#' `lambdatwo`) to the congruence–bias parameterisation (the midpoint `M` of the two
#' expected counts, the congruence `f`, and the source bias `delta`), and,
#' optionally, converts native dispersions (`shapes`, `shapex`, `shapexone`,
#' and `shapextwo`) to overdispersions (`kappas`, `kappax`, `kappaxone`, and
#' `kappaxtwo`). [mfd_to_dpars()] converts in the other direction.
#'
#' @param mu Rate of the shared component.
#' @param lambdaone,lambdatwo Rates of the source-1-exclusive and
#'   source-2-exclusive components.
#' @param shapes,shapex Optional NB2 (inverse) dispersion parameters.
#'   Supplying `shapes` adds `kappas = 1/sqrt(shapes)` to the returned list, so
#'   `shapes = Inf` gives `kappas = 0`. A `binegbin` fit made before 0.10.0
#'   has one dispersion, `shapex`, for both exclusive components; supplying
#'   `shapex` adds `kappax`. For [bipois()] and
#'   [bipois_partialobs()], which have no dispersion parameters, leave every
#'   dispersion argument unset.
#' @param shapexone,shapextwo Optional NB2 (inverse) dispersion parameters of
#'   the two exclusive components. Supplying them adds `kappaxone` and
#'   `kappaxtwo` to the returned list. Mutually exclusive with `shapex`.
#'
#' @details
#' Each argument may be a vector, such as the posterior draws of the rates.
#' `mu`, `lambdaone`, and `lambdatwo` are recycled to the length of the
#' longest of the three, and any dispersions to the same length; each element
#' of the returned list is a vector of that length.
#'
#' Rates of exactly zero, which posterior draws do not produce, give values
#' that are not finite: `delta` is `NA` when both exclusive rates are zero,
#' and `Inf` or `-Inf` when one of them is, and `f` is `NA` when `M` is zero.
#'
#' @return A named list of `M`, `f`, `delta`, and `beta = tanh(delta)`, the
#'   source bias on a scale bounded by -1 and 1, with `kappas`, `kappax`,
#'   `kappaxone`, and `kappaxtwo` added when the corresponding dispersions are
#'   supplied.
#'
#' @examples
#' dpars_to_mfd(mu = 8.04, lambdaone = 4.75, lambdatwo = 3.17)
#'
#' # Round trip
#' d <- mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#' dpars_to_mfd(d$mu, d$lambdaone, d$lambdatwo)[c("M", "f", "delta")]
#'
#' # Perfect congruence: delta cannot be recovered and is returned as NA
#' dpars_to_mfd(mu = 12, lambdaone = 0, lambdatwo = 0)$delta
#'
#' @seealso [mfd_to_dpars()], [bipois()], [bipois_partialobs()],
#'   [binegbin()], [binegbin_partialobs()]
#' @export
dpars_to_mfd <- function(mu, lambdaone, lambdatwo,
                         shapes = NULL, shapex = NULL,
                         shapexone = NULL, shapextwo = NULL) {
  n <- max(length(mu), length(lambdaone), length(lambdatwo))
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  if (any(c(mu, lambdaone, lambdatwo) < 0, na.rm = TRUE)) {
    stop("Rates must be non-negative.", call. = FALSE)
  }
  # Mirrors the guard in mfd_to_dpars(): one spelling of the
  # exclusive-component dispersion per call, chosen by how old the fit is.
  if (!is.null(shapex) && (!is.null(shapexone) || !is.null(shapextwo))) {
    stop("Supply either `shapex` (one dispersion for both exclusive ",
         "components, from a binegbin fit made before 0.10.0) or ",
         "`shapexone`/`shapextwo` (a dispersion for each exclusive ",
         "component), not both.", call. = FALSE)
  }

  exclusive_sum <- lambdaone + lambdatwo
  M <- mu + exclusive_sum / 2

  # f undefined when there is nothing at all.
  f <- ifelse(M > 0, mu / M, NA_real_)

  # delta undefined when there are no exclusive events to be biased (f == 1).
  # One rate zero and the other positive is the legitimate +/-Inf limit, which
  # log() produces directly.
  delta <- ifelse(exclusive_sum > 0, 0.5 * log(lambdaone / lambdatwo), NA_real_)
  beta  <- ifelse(exclusive_sum > 0, (lambdaone - lambdatwo) / exclusive_sum, NA_real_)

  out <- list(M = M, f = f, delta = delta, beta = beta)

  if (!is.null(shapes))    out$kappas    <- .shape_to_kappa(rep_len(shapes, n))
  if (!is.null(shapex))    out$kappax    <- .shape_to_kappa(rep_len(shapex, n))
  if (!is.null(shapexone)) out$kappaxone <- .shape_to_kappa(rep_len(shapexone, n))
  if (!is.null(shapextwo)) out$kappaxtwo <- .shape_to_kappa(rep_len(shapextwo, n))
  out
}

# Internal conversions between inverse dispersion and overdispersion. Note the
# direction reversal: kappa increases with overdispersion, shape decreases.
# kappa = 0 <-> shape = Inf is the Poisson limit and is handled exactly rather
# than by division blowing up.
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
