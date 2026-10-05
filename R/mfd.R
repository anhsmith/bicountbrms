# ==========================================================================
# (M, f, delta) <-> native dpar coordinates
#
# THE RATE HALF SERVES BOTH FAMILIES, AND BOTH CONSTRUCTORS OF EACH. mu,
# lambdaone and lambdatwo are the same three rates in bipois(),
# bipois_partialobs(), binegbin() and binegbin_partialobs(), so the
# (M, f, delta) map below applies to every one of them unchanged. Only the
# DISPERSION half differs:
#
#   bipois(), bipois_partialobs()   no dispersion -- supply no kappa
#   binegbin(), binegbin_partialobs()
#                                   one per margin,
#                                   kappaxone/kappaxtwo -> shapexone/shapextwo
#
# THE TWO DIRECTIONS SPEAK DIFFERENT VOCABULARIES, DELIBERATELY.
# binegbin_mfd_to_dpars() emits only names a SHIPPING family accepts, because
# its output is meant to be pasted into a brm() call or a simulation. Since
# 0.10.0 no family has a `shapex` dpar, so `kappax` -- the shorthand for "one
# overdispersion for both exclusive components" -- writes shapexone and
# shapextwo at the same value rather than a single `shapex`.
#
# binegbin_dpars_to_mfd() still ACCEPTS `shapex`, because its input is a
# stored fit and `shapex` is the genuine dpar name on any binegbin fit made
# before 0.10.0. This is the same principle .get_dpar_any() applies in
# utils.R: read the old spellings, write the current ones.
#
# binegbin()/bipois() are parameterised by three rates -- mu (shared),
# lambdaone and lambdatwo (the rates of the two source-exclusive components) --
# because the trivariate-reduction likelihood is written in terms of these three
# rates. mu, lambdaone and lambdatwo are correlated in use: increasing the size
# of the catch moves all three at once, so none of the three rates alone answers
# "what is the midpoint of the two expected counts", "what proportion of the
# midpoint did both sources record", or "which source recorded more of the
# unshared events".
#
# The (M, f, delta) coordinates separate exactly those three questions:
#
#   M     = mu + (lambdaone + lambdatwo)/2  midpoint of the two expected
#                                           counts, (E[y1] + E[y2])/2
#   f     = mu / M                          congruence: the mean of the shared
#                                           component as a proportion of M
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
# as the fractional imbalance between the two exclusive rates, and
# beta = +/-1 is the limit where one source never records an unshared event.
#
# The map is a bijection on the interior. Because the average of the two
# exclusive rates is M(1-f) for ANY delta, M stays pinned to the midpoint
# whatever the bias -- M is a midpoint of what the two sources REPORT, not a
# property of the underlying process.
#
# INVERSE DISPERSION AND OVERDISPERSION. The dpars shapes/shapexone/shapextwo
# are NB2 phi (Stan neg_binomial_2, R dnbinom size): variance = m + m^2/phi, so
# LARGER phi means LESS overdispersion and Poisson is the phi -> Inf limit. The
# overdispersion kappa, used alongside (M, f, delta), inverts that into a dial
# that increases with overdispersion and reaches Poisson at a finite zero:
#
#   shapes = 1/kappas^2      kappas = 1/sqrt(shapes)
#
# Note the direction reversal: raising kappa LOWERS shape.
#
# The two converters in this file are pure coordinate transforms -- they fit
# nothing. To fit in (M, f, delta) coordinates, supply them through a
# non-linear formula; see the examples on binegbin_mfd_to_dpars(). Use the
# converters to set up simulations from interpretable values, and to read
# fitted dpars back into the interpretable coordinates.

# --------------------------------------------------------------------------

#' Convert (M, f, delta) coordinates to native binegbin/bipois dpars
#'
#' @description
#' Maps the interpretable coordinates -- the midpoint `M` of the two expected
#' counts, congruence `f`, and source bias `delta` -- onto the rate dpars every
#' family in this package takes (`mu`, `lambdaone`, `lambdatwo`), optionally
#' converting overdispersions to the `shapes`/`shapexone`/`shapextwo`
#' dpars.
#'
#' The three rates are common to [bipois()], [bipois_partialobs()],
#' [binegbin()] and [binegbin_partialobs()], so this direction serves all four
#' constructors. The four constructors differ in their dispersions: the
#' Poisson families have none, and both negative-binomial constructors have one
#' per margin.
#'
#' Everything this returns is named for a dpar a shipping family accepts, so
#' the output can go straight into a `brm()` call or a simulation. That is why
#' `kappax` writes `shapexone` and `shapextwo` at a common value rather than a
#' single `shapex`: no family has had a `shapex` dpar since 0.10.0.
#'
#' [binegbin_dpars_to_mfd()] is the inverse. It still accepts `shapex`, because
#' its input is a stored fit and pre-0.10.0 `binegbin` fits declare that name.
#'
#' @param M Midpoint of the two expected counts:
#'   `mu + (lambdaone + lambdatwo)/2`, equal to `(E[y1] + E[y2])/2`.
#'   Non-negative.
#' @param f Congruence, the mean of the shared component as a proportion of
#'   `M`: `mu / M`. In `[0, 1]`. `f = 1` means perfect congruence (both
#'   source-exclusive components vanish); `f = 0` means no shared component at
#'   all.
#' @param delta Source bias on the log-ratio scale,
#'   `0.5 * log(lambdaone/lambdatwo)`. `0` is unbiased. `+/-Inf` is permitted
#'   and gives the limit where one exclusive rate is zero.
#' @param kappas,kappax Optional overdispersions, `kappa = 1/sqrt(phi)`. `0` is
#'   the Poisson limit. `kappas` is the overdispersion of the shared component,
#'   and the returned list gains `shapes` (`= 1/kappa^2`, so `kappa = 0` gives
#'   `Inf`). `kappax` is the shorthand for a single overdispersion governing
#'   *both* margins: supply it and the returned list gains `shapexone` and
#'   `shapextwo` at that common value, which is the symmetric model [binegbin()]
#'   reaches by tying the two with `nlf()`. Omit both for [bipois()] and
#'   [bipois_partialobs()], which have no dispersion parameters.
#' @param kappaxone,kappaxtwo Optional per-margin overdispersions of the
#'   exclusive components, for the general case in which the two margins are
#'   free to differ. If supplied, the returned list gains
#'   `shapexone`/`shapextwo`. Mutually exclusive with `kappax`, which writes the
#'   same two slots.
#'
#' @details
#' Arguments are recycled to a common length, so this vectorises over posterior
#' draws.
#'
#' **Boundary behaviour.** At `f = 1` both exclusive rates are exactly `0`
#' regardless of `delta` -- the bias becomes unidentifiable, which
#' [binegbin_dpars_to_mfd()] reports back as `NA`. This direction is always
#' well defined; only the inverse degenerates.
#'
#' **The two spellings of the overdispersion of the exclusive components.** Both
#' negative-binomial constructors take the pair `shapexone`/`shapextwo`, so this
#' function returns `shapexone` and `shapextwo` and never `shapex`. Supply
#' `kappaxone`/`kappaxtwo` to give the two margins different values, or
#' `kappax` to give both margins a single value -- the symmetric model, term
#' for term the pre-0.8.0 likelihood, which a fit reaches by routing both dpars
#' through one non-linear parameter. The two spellings are mutually exclusive
#' in a single call because they write the same two slots.
#'
#' Before 0.10.0, `kappax` returned a dpar named `shapex` and
#' `kappaxone`/`kappaxtwo` returned the pair, because two different families
#' wanted two different things. There is now one family and one dpar set, so
#' both spellings produce it.
#'
#' **These coordinates under the Poisson families.** [bipois()] and
#' [bipois_partialobs()] take the same three rates and no dispersion, so
#' call this with `M`, `f` and `delta` alone and pass the result straight
#' through. There is no `kappa` to supply: the Poisson case is not an
#' dispersion set to a particular value but the absence of the parameter, which
#' is precisely why fitting it wants its own family rather than [binegbin()]
#' with `kappa` driven to `0`.
#'
#' @return A named list of `mu`, `lambdaone`, `lambdatwo`, plus `shapes` when
#'   `kappas` is supplied, and `shapexone`/`shapextwo` when either `kappax` or
#'   `kappaxone`/`kappaxtwo` are supplied. Every name is a dpar of a shipping
#'   family.
#'
#' @examples
#' # A moderately congruent pair, source 1 recording more of the unshared
#' # events
#' binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#'
#' # Perfect congruence: both exclusive components vanish
#' binegbin_mfd_to_dpars(M = 12, f = 1, delta = 0.5)
#'
#' # kappax: one overdispersion for both exclusive components
#' binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4, kappax = 0.9)
#'
#' # kappaxone, kappaxtwo: one overdispersion per exclusive component
#' binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4,
#'                       kappaxone = 0.9, kappaxtwo = 0.3)
#'
#' # To FIT in these coordinates, pass them through a non-linear formula
#' # (every dpar is log-linked, so the link supplies the exp()):
#' #   bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#' #     nlf(mu        ~ logM + log_inv_logit(logitf)) +
#' #     nlf(lambdaone ~ log(2) + logM + log_inv_logit(-logitf) +
#' #                     log_inv_logit(2 * delta)) +
#' #     nlf(lambdatwo ~ log(2) + logM + log_inv_logit(-logitf) +
#' #                     log_inv_logit(-2 * delta)) +
#' #     lf(logM ~ 1, logitf ~ 1, delta ~ 1)
#' # where logM = log M, logitf = logit f and delta = artanh(beta).
#'
#' @seealso [binegbin_dpars_to_mfd()], [bipois()], [bipois_partialobs()],
#'   [binegbin()], [binegbin_partialobs()]
#' @export
binegbin_mfd_to_dpars <- function(M, f, delta = 0, kappas = NULL, kappax = NULL,
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
  # `shapex` here would name a dpar no shipping family accepts, and the output
  # of this function is meant to be usable in a brm() call unaltered.
  if (!is.null(kappax)) {
    # Assigned separately, not chained: `a <- b <- value` evaluates
    # right-to-left, so the chained form would create shapextwo first and the
    # returned list would read shapextwo, shapexone. Dpar order is load-bearing
    # in this package, and a reader may infer the wrong signature from a
    # printed list in the wrong order.
    kx <- .kappa_to_shape(rep_len(kappax, n))
    out$shapexone <- kx
    out$shapextwo <- kx
  }
  if (!is.null(kappaxone)) out$shapexone <- .kappa_to_shape(rep_len(kappaxone, n))
  if (!is.null(kappaxtwo)) out$shapextwo <- .kappa_to_shape(rep_len(kappaxtwo, n))
  out
}

#' Convert native binegbin/bipois dpars to (M, f, delta) coordinates
#'
#' @description
#' Inverse of [binegbin_mfd_to_dpars()]. Reads the rate dpars `mu`,
#' `lambdaone`, `lambdatwo` back into the interpretable coordinates: the
#' midpoint `M` of the two expected counts, congruence `f`, and source bias
#' `delta`, optionally converting NB2 inverse dispersions to overdispersions.
#'
#' As with the forward direction, the three rates are common to [bipois()],
#' [bipois_partialobs()], [binegbin()] and [binegbin_partialobs()], so this
#' serves all four constructors; only the dispersion arguments differ.
#'
#' This direction reads a *stored fit*, so it accepts `shapex` -- the single
#' dispersion of both exclusive components, declared by every
#' `binegbin` fit made before 0.10.0 -- as well as the current
#' `shapexone`/`shapextwo`. The forward direction writes only current names;
#' see its documentation for why the two directions differ.
#'
#' @param mu Shared-component rate.
#' @param lambdaone,lambdatwo The rates of the two source-exclusive components.
#' @param shapes,shapex Optional NB2 inverse dispersions. If supplied, the
#'   returned list gains `kappas`/`kappax` (`= 1/sqrt(shape)`, so `shape = Inf`
#'   gives `0`). `shapex` is the single dispersion of both exclusive
#'   components, declared by a `binegbin` fit made before 0.10.0;
#'   no shipping family takes it now. A fit of [bipois()] or
#'   [bipois_partialobs()] has no dispersion to pass: read its three
#'   rates alone.
#' @param shapexone,shapextwo Optional NB2 inverse dispersions of the exclusive
#'   components, as taken by both negative-binomial constructors since 0.10.0.
#'   If supplied, the returned list gains `kappaxone`/`kappaxtwo`. Mutually
#'   exclusive with `shapex`, which is the same quantity under its older name.
#'
#' @details
#' Arguments are recycled to a common length, so this vectorises over posterior
#' draws -- e.g. to convert a whole posterior into interpretable coordinates.
#'
#' **Boundary behaviour**, which the forward direction does not have:
#'
#' * `lambdaone == lambdatwo == 0` (perfect congruence, `f = 1`): the bias is
#'   genuinely unidentified -- there are no exclusive events to be biased -- and
#'   `delta` is returned as `NA`, not `0`. Zero would assert an unbiased source,
#'   which the data at that point cannot support.
#' * `M == 0` (nothing anywhere): `f` is undefined and returned as `NA`.
#' * Exactly one exclusive rate `0`: `delta` is `+/-Inf`, the well-defined limit
#'   where one source never records an unshared event.
#'
#' Because of the first case, round-tripping is exact everywhere except at
#' `f = 1`, where `delta` cannot be recovered. Hold one coordinate system as the
#' source of truth rather than repeatedly converting back and forth.
#'
#' @return A named list of `M`, `f`, `delta`, plus `beta` (the bounded bias
#'   `tanh(delta)`), `kappas`/`kappax` when `shapes`/`shapex` are supplied, and
#'   `kappaxone`/`kappaxtwo` when `shapexone`/`shapextwo` are supplied.
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
#' @seealso [binegbin_mfd_to_dpars()], [bipois()], [bipois_partialobs()],
#'   [binegbin()], [binegbin_partialobs()]
#' @export
binegbin_dpars_to_mfd <- function(mu, lambdaone, lambdatwo,
                                  shapes = NULL, shapex = NULL,
                                  shapexone = NULL, shapextwo = NULL) {
  n <- max(length(mu), length(lambdaone), length(lambdatwo))
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  if (any(c(mu, lambdaone, lambdatwo) < 0, na.rm = TRUE)) {
    stop("Rates must be non-negative.", call. = FALSE)
  }
  # Mirrors the guard in binegbin_mfd_to_dpars(): one spelling of the
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
