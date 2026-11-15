#' Random generation of paired counts
#'
#' @description
#' `rbicount()` draws pairs of counts `(y1, y2)` from the generating model of
#' the `bipois` and `binegbin` families, and returns the three latent
#' components of each pair with the pair. The arguments are the native
#' parameters of the families, in the order and under the names that
#' [mfd_to_dpars()] returns them.
#'
#' @param n Number of pairs to draw. If `length(n) > 1`, the length of `n` is
#'   taken to be the number required, as in [stats::rnbinom()].
#' @param mu Mean of the shared component `Ns`.
#' @param lambdaone,lambdatwo Means of the source-exclusive components `N1`
#'   and `N2`.
#' @param shapes,shapexone,shapextwo NB2 (inverse) dispersion parameters of
#'   `Ns`, `N1`, and `N2`. Supply all three for negative-binomial components
#'   (the `binegbin` family), or none for Poisson components (the `bipois`
#'   family).
#'
#' @details
#' Both families construct the pair by trivariate reduction. Three independent
#' latent counts are drawn: the shared component `Ns`, recorded by both
#' sources, with mean `mu`; and the source-exclusive components `N1` and `N2`,
#' each recorded by one source alone, with means `lambdaone` and `lambdatwo`.
#' The recorded counts are
#'
#' \deqn{y_1 = N_s + N_1, \qquad y_2 = N_s + N_2,}{y1 = Ns + N1,
#'   y2 = Ns + N2,}
#'
#' so that \eqn{E[y_1] = \mu + \lambda_1}{E[y1] = mu + lambdaone},
#' \eqn{E[y_2] = \mu + \lambda_2}{E[y2] = mu + lambdatwo}, and
#' \eqn{\mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_s)}{Cov(y1, y2) = Var(Ns)}.
#'
#' With the three dispersions supplied, each component is negative binomial
#' in the NB2 parameterisation, drawn by `rnbinom(size = phi, mu = m)`. A
#' component with mean \eqn{m}{m} and dispersion \eqn{\phi}{phi} has variance
#' \eqn{m + m^2/\phi}{m + m^2/phi}, so a small \eqn{\phi}{phi} gives a heavy
#' tail and \eqn{\phi = \infty}{phi = Inf} gives the Poisson limit. The
#' dispersion of `Ns` is `shapes`, of `N1` is `shapexone`, and of `N2` is
#' `shapextwo`. With all three dispersions `NULL`, each component is Poisson,
#' with variance equal to its mean. Supplying some of the dispersions and not
#' others stops with an error, because neither family has a mixture of Poisson
#' and negative-binomial components.
#'
#' Each parameter is recycled to length `n`, as by [stats::rnbinom()], so a
#' parameter may be a single value or a vector with one value per pair. The
#' list returned by [mfd_to_dpars()] has elements named `mu`, `lambdaone`,
#' `lambdatwo`, and, when overdispersions are supplied to it, `shapes`,
#' `shapexone`, and `shapextwo`, so `do.call(rbicount, c(list(n = n), pars))`
#' draws pairs at the native parameters `pars` (see Examples). The overdispersion
#' \eqn{\kappa = 0}{kappa = 0} gives `shapes = Inf` in [mfd_to_dpars()], and
#' `rbicount()` then draws a Poisson component.
#'
#' The components are drawn by the code that [posterior_predict_y1_and_y2()]
#' uses for each posterior draw, so a pair simulated by `rbicount()` and a
#' posterior predictive draw at the same parameter values have the same
#' distribution.
#'
#' To simulate a partially paired design, draw an observation indicator, such as
#' `y1_obs <- rbinom(n, 1, 0.2)`, and set `y1` to `-1` on the rows
#' where `y1_obs` is 0 before fitting with [bipois_partialobs()] or
#' [binegbin_partialobs()]. The columns `Ns` and `N1` keep the values from
#' which the unrecorded `y1` would have been formed.
#'
#' @return A data frame with `n` rows and five integer columns: `y1` and `y2`,
#'   the two recorded counts, followed by `Ns`, `N1`, and `N2`, the shared
#'   component and the two exclusive components, with `y1 = Ns + N1` and
#'   `y2 = Ns + N2` on every row.
#'
#' @seealso [mfd_to_dpars()] for the native parameters at given values of the
#'   midpoint, congruence, and source bias; [posterior_predict_y1_and_y2()]
#'   for draws of the pair from a fitted model.
#'
#' @examples
#' set.seed(1)
#'
#' # Poisson components, as in bipois()
#' rbicount(5, mu = 4, lambdaone = 1, lambdatwo = 2)
#'
#' # Negative-binomial components, as in binegbin(). The covariance of the
#' # two counts is Var(Ns) = 4 + 4^2/2 = 12.
#' d <- rbicount(
#'   10000,
#'   mu        = 4,
#'   lambdaone = 1,
#'   lambdatwo = 2,
#'   shapes    = 2,
#'   shapexone = 4,
#'   shapextwo = 4
#' )
#' cov(d$y1, d$y2)
#'
#' # Parameters from the midpoint, congruence, and source bias
#' pars <- mfd_to_dpars(
#'   M      = 5,
#'   f      = 0.7,
#'   delta  = atanh(0.3),
#'   kappas = 0.5,
#'   kappax = 0.5
#' )
#' d <- do.call(rbicount, c(list(n = 100), pars))
#'
#' # A partially paired design. y2 is recorded on every unit, and y1 on about
#' # 20% of units.
#' d$y1_obs  <- rbinom(nrow(d), 1, 0.2)
#' d$y1_true <- d$y1
#' d$y1[d$y1_obs == 0] <- -1L
#' head(d)
#'
#' # Thirty units in each of three treatments that differ in congruence
#' if (requireNamespace("dplyr", quietly = TRUE)) {
#'   units <- data.frame(
#'     treatment = rep(c("a", "b", "c"), each = 30)
#'   )
#'   units |>
#'     dplyr::mutate(
#'       f = c(a = 0.5, b = 0.7, c = 0.9)[treatment],
#'       as.data.frame(
#'         mfd_to_dpars(
#'           M      = 5,
#'           f      = f,
#'           delta  = 0,
#'           kappas = 0.5,
#'           kappax = 0.5
#'         )
#'       ),
#'       rbicount(
#'         dplyr::n(),
#'         mu,
#'         lambdaone,
#'         lambdatwo,
#'         shapes,
#'         shapexone,
#'         shapextwo
#'       )
#'     ) |>
#'     head()
#' }
#'
#' @export
rbicount <- function(n,
                     mu,
                     lambdaone,
                     lambdatwo,
                     shapes    = NULL,
                     shapexone = NULL,
                     shapextwo = NULL) {
  if (length(n) > 1L) {
    n <- length(n)
  }
  is_null <- c(
    shapes    = is.null(shapes),
    shapexone = is.null(shapexone),
    shapextwo = is.null(shapextwo)
  )
  if (any(is_null) && !all(is_null)) {
    stop(
      "Supply all three dispersions (`shapes`, `shapexone`, and `shapextwo`) ",
      "for negative-binomial components, or none of them for Poisson ",
      "components. Missing: ",
      paste0("`", names(is_null)[is_null], "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  comp <- .draw_components(
    n,
    mu,
    lambdaone,
    lambdatwo,
    shapes    = shapes,
    shapexone = shapexone,
    shapextwo = shapextwo
  )
  # rnbinom() returns doubles and rpois() integers. Integer columns under both
  # constructors let a placeholder such as -1L be assigned to y1 without a
  # change of type, which dplyr::if_else() would reject.
  comp <- lapply(comp, as.integer)
  data.frame(
    y1 = comp$Ns + comp$N1,
    y2 = comp$Ns + comp$N2,
    Ns = comp$Ns,
    N1 = comp$N1,
    N2 = comp$N2
  )
}
