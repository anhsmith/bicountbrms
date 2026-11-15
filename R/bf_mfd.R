# ==========================================================================
# bf_mfd(): the brms non-linear formula in (M, f, delta, kappa) coordinates
#
#   mu        = M f
#   lambdaone = M (1 - f)(1 + tanh delta) = 2 M (1 - f) inv_logit(2 delta)
#   lambdatwo = M (1 - f)(1 - tanh delta) = 2 M (1 - f) inv_logit(-2 delta)
#   shape     = 1 / kappa^2
#
# written on the log scale of the native dpars, which all four constructors
# log-link. The rate lines are the same identities as mfd_to_dpars() in
# R/mfd.R; test-bf-mfd.R checks the two against each other. The function only
# builds a formula: it adds no family, Stan function or post-processing
# method, because the fitted model is an ordinary bipois or binegbin fit whose
# dpars are non-linear functions of the coordinates.
#
# The rate lines use the inv_logit form of the source-bias factor. The term
# log_inv_logit(-2 * delta) is exact for any delta, whereas
# log1m(tanh(delta)) loses precision as tanh(delta) approaches 1 and returns
# -Inf from delta = 20. Under a weak prior on delta, the sampler can reach
# such values when one source records almost no exclusive events.

.MFD_RATE_LINES <- c(
  mu        = "mu ~ logM + log_inv_logit(logitf)",
  lambdaone = "lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)",
  lambdatwo = "lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)"
)

.MFD_COORDS <- c("logM", "logitf", "delta",
                 "kappas", "kappax", "kappaxone", "kappaxtwo")

# The three valid sets of coordinate formulas, and the nlf() lines each adds to
# the rate lines.
.MFD_SETS <- list(
  separate = list(
    coords = c("logM", "logitf", "delta", "kappas", "kappaxone", "kappaxtwo"),
    lines  = c(shapes    = "shapes ~ -2 * log(kappas)",
               shapexone = "shapexone ~ -2 * log(kappaxone)",
               shapextwo = "shapextwo ~ -2 * log(kappaxtwo)")
  ),
  tied = list(
    coords = c("logM", "logitf", "delta", "kappas", "kappax"),
    lines  = c(shapes    = "shapes ~ -2 * log(kappas)",
               shapexone = "shapexone ~ -2 * log(kappax)",
               shapextwo = "shapextwo ~ -2 * log(kappax)")
  ),
  poisson = list(
    coords = c("logM", "logitf", "delta"),
    lines  = character(0)
  )
)

#' Build a brms formula in the congruence–bias parameterisation
#'
#' @description
#' Builds the brms non-linear formula for the congruence–bias parameterisation,
#' in which the native parameters of a [bipois()] or [binegbin()] model are
#' functions of the congruence–bias coordinates: the midpoint `M` of the two
#' expected counts, the congruence `f`, the source bias (through `delta`, its
#' inverse hyperbolic tangent), and, for `binegbin`, the overdispersions
#' `kappa`. Each coordinate then takes its own linear predictor. A predictor
#' such as site in the formula for `logM` alone multiplies all three rates by a
#' common factor, so the congruence and the source bias are the same at every
#' site. In the native parameterisation, the same model needs the predictor in
#' all three rates, with the coefficients constrained to be equal.
#'
#' @param formula The response part of the model, written without `~`:
#'   `y1 | vint(y2)` for [bipois()] and [binegbin()], and
#'   `y1 | vint(y2, y1_obs)` for [bipois_partialobs()] and
#'   [binegbin_partialobs()]. The response is passed to the brms formula
#'   unchanged.
#' @param ... Two-sided formulas, one for each coordinate, in any order. The
#'   left-hand side names the coordinate and the right-hand side gives its
#'   linear predictor, as in `logM ~ 0 + site` or `delta ~ 1`. Every
#'   coordinate of the chosen set is required, including those with an
#'   intercept only.
#'
#' @details
#' The coordinates are related to the native parameters by
#'
#' \deqn{\mu = M f, \quad \lambda_1 = M (1 - f)(1 + \beta), \quad
#'   \lambda_2 = M (1 - f)(1 - \beta), \quad \beta = \tanh \delta,}
#'
#' where \eqn{\mu} is `mu`, the rate of the shared component, and
#' \eqn{\lambda_1} and \eqn{\lambda_2} are `lambdaone` and `lambdatwo`, the
#' rates of the source-exclusive components. Each NB2 (inverse) dispersion
#' parameter \eqn{\phi} is \eqn{1/\kappa^2}, where \eqn{\kappa} is the
#' corresponding overdispersion. The non-linear parameters are `logM`
#' (\eqn{\log M}), `logitf` (\eqn{\operatorname{logit} f}), `delta`
#' (\eqn{\delta}), and the overdispersions. `logM`, `logitf`, and `delta` are
#' unbounded, so their linear predictors need no link. The overdispersions
#' must be positive (see Priors). [mfd_to_dpars()]
#' computes the native parameters from the same identities.
#'
#' **Sets of coordinates.** The set of coordinate formulas supplied determines
#' the dispersion structure. `logM`, `logitf`, and `delta` are always
#' required. The three valid sets are:
#'
#' * `logM`, `logitf`, `delta`, `kappas`, `kappaxone`, and `kappaxtwo`: a
#'   separate overdispersion for the shared component and for each exclusive
#'   component, for [binegbin()] and [binegbin_partialobs()].
#' * `logM`, `logitf`, `delta`, `kappas`, and `kappax`: one overdispersion,
#'   `kappax`, for both exclusive components, so that `shapexone` equals
#'   `shapextwo`, for [binegbin()] and [binegbin_partialobs()].
#' * `logM`, `logitf`, and `delta`: no overdispersion, for [bipois()] and
#'   [bipois_partialobs()].
#'
#' Any other set stops with an error. The family passed to [brms::brm()] must
#' match the set: a set with overdispersions written for a `bipois` model
#' declares dpars that the family does not have.
#'
#' **The generated formula.** For the second set with a site effect on `logM`,
#'
#' ```
#' bf_mfd(
#'   y1 | vint(y2),
#'   logM   ~ 0 + site,
#'   logitf ~ 1,
#'   delta  ~ 1,
#'   kappas ~ 1,
#'   kappax ~ 1
#' )
#' ```
#'
#' returns the same formula as
#'
#' ```
#' bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#'   nlf(mu        ~ logM + log_inv_logit(logitf)) +
#'   nlf(lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)) +
#'   nlf(lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta)) +
#'   nlf(shapes    ~ -2 * log(kappas)) +
#'   nlf(shapexone ~ -2 * log(kappax)) +
#'   nlf(shapextwo ~ -2 * log(kappax)) +
#'   lf(logM ~ 0 + site, logitf ~ 1, delta ~ 1, kappas ~ 1, kappax ~ 1)
#' ```
#'
#' The first set replaces the two `kappax` lines with
#' `shapexone ~ -2 * log(kappaxone)` and `shapextwo ~ -2 * log(kappaxtwo)`;
#' the third set omits all three dispersion lines.
#'
#' **Priors.** Every coordinate is a non-linear parameter, to which brms gives
#' a flat, improper prior unless one is set. We recommend a prior on each,
#' written with `class = "b"` and `nlpar` set to the name of the coordinate.
#' The expression \eqn{-2 \log \kappa} is undefined for
#' \eqn{\kappa \le 0}, so a prior on an overdispersion with an intercept-only
#' formula needs `lb = 0`, as in
#' `prior(normal(0, 1), class = "b", nlpar = "kappas", lb = 0)`. The worked
#' example in `vignette("bicountbrms")` gives priors
#' for every coordinate.
#'
#' `bf_mfd()` adds no family, Stan function, or post-processing method. The
#' fitted model is a `bipois` or `binegbin` fit whose native parameters are
#' non-linear functions of the coordinates, and `posterior_epred()`,
#' `posterior_predict()`, and `log_lik()` dispatch to the family methods as for
#' any other fit.
#'
#' @return A `brmsformula` object, as returned by [brms::bf()].
#'
#' @examples
#' # A site effect on the midpoint alone, with one overdispersion for both
#' # exclusive components
#' bf_mfd(
#'   y1 | vint(y2),
#'   logM   ~ 0 + site,
#'   logitf ~ 1,
#'   delta  ~ 1,
#'   kappas ~ 1,
#'   kappax ~ 1
#' )
#'
#' # A partially paired design with Poisson components
#' bf_mfd(
#'   y1 | vint(y2, y1_obs),
#'   logM   ~ 0 + site,
#'   logitf ~ 1,
#'   delta  ~ 1
#' )
#'
#' @seealso [mfd_to_dpars()], [dpars_to_mfd()], [bipois()],
#'   [bipois_partialobs()], [binegbin()], [binegbin_partialobs()]
#' @export
bf_mfd <- function(formula, ...) {
  if (missing(formula)) {
    stop("Supply the response as the first argument, as in ",
         "`bf_mfd(y1 | vint(y2), logM ~ 1, ...)`.", call. = FALSE)
  }
  env  <- parent.frame()
  resp <- substitute(formula)
  # A formula in the first position is either the response written with `~`
  # or a coordinate formula given before the response. Neither can be passed
  # through, and the second would otherwise reach brms as a response named
  # `logM`.
  if (is.call(resp) && identical(resp[[1]], as.name("~"))) {
    stop("The first argument of `bf_mfd()` is the response, written without ",
         "`~`, as in `y1 | vint(y2)` or `y1 | vint(y2, y1_obs)`. ",
         "The coordinate formulas follow it.", call. = FALSE)
  }

  coord_forms <- list(...)
  is_twosided <- vapply(
    coord_forms,
    function(x) inherits(x, "formula") && length(x) == 3L,
    logical(1)
  )
  if (!all(is_twosided)) {
    stop("Every argument after the response must be a two-sided formula ",
         "naming a coordinate, as in `logM ~ 0 + site` or `delta ~ 1`.",
         call. = FALSE)
  }

  lhs <- vapply(coord_forms, function(x) deparse(x[[2]]), character(1))
  unknown <- setdiff(lhs, .MFD_COORDS)
  if (length(unknown) > 0) {
    stop("Unknown coordinate ", paste0("`", unknown, "`", collapse = ", "),
         " on the left-hand side of a formula. The coordinates are ",
         paste0("`", .MFD_COORDS, "`", collapse = ", "), ".", call. = FALSE)
  }
  dup <- unique(lhs[duplicated(lhs)])
  if (length(dup) > 0) {
    stop("More than one formula for ",
         paste0("`", dup, "`", collapse = ", "), ".", call. = FALSE)
  }

  set <- NULL
  for (s in .MFD_SETS) {
    if (setequal(lhs, s$coords)) set <- s
  }
  if (is.null(set)) {
    stop("The coordinate formulas supplied (",
         paste0("`", lhs, "`", collapse = ", "),
         ") are not one of the three valid sets:\n",
         "  * logM, logitf, delta, kappas, kappaxone, kappaxtwo ",
         "(binegbin, separate overdispersions of the exclusive components);\n",
         "  * logM, logitf, delta, kappas, kappax ",
         "(binegbin, one overdispersion for both exclusive components);\n",
         "  * logM, logitf, delta (bipois).", call. = FALSE)
  }

  # Canonical order, so that the generated formula, and the Stan code, do not
  # depend on the order of the arguments.
  coord_forms <- coord_forms[match(set$coords, lhs)]

  out <- brms::bf(eval(call("~", resp, 1), env), nl = TRUE)
  for (line in c(.MFD_RATE_LINES, set$lines)) {
    out <- out + brms::nlf(stats::as.formula(line, env = env))
  }
  out + do.call(brms::lf, coord_forms)
}
