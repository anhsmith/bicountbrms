# ==========================================================================
# Prediction of the pair, and imputation of the first count
# ==========================================================================
#
# Both families generate a pair of counts by trivariate reduction,
#
#   y1 = N_shared + N1,   y2 = N_shared + N2,
#
# with N_shared, N1, and N2 independent given their means mu, lambdaone, and
# lambdatwo. brms supplies y2 to the family as data through vint(), so a
# post-processing method could either condition on the recorded y2 or ignore
# the recorded y2.
#
# The brms methods ignore the recorded y2. posterior_predict_<family>() draws
# y1 from its marginal distribution given the parameters, and
# posterior_epred_<family>() returns the marginal mean mu + lambdaone, so
# that pp_check(), fitted(), conditional_effects(), and bayes_R2() have the
# meaning they have for any other brms family. This file adds four functions,
# in the brms convention of one function for draws and one for expectations:
#
#   posterior_predict_y1_and_y2()    draws of (y1, y2) from the joint model
#   posterior_epred_y1_and_y2()      (mu + lambdaone, mu + lambdatwo)
#   posterior_predict_y1_given_y2()  draws from p(y1 | y2, theta)
#   posterior_epred_y1_given_y2()    E[y1 | y2] = E[N_shared | y2] + lambdaone
#
# brms finds a custom-family method by pasting posterior_predict_,
# posterior_epred_, or log_lik_ onto the family name and looking the result up
# in the environment of the family. The family names are bipois and binegbin,
# so these four names cannot be mistaken for methods unless a family is ever
# named `y1_and_y2` or `y1_given_y2`.
#
# The prep-level workers (.impute_draw_<family>(), .impute_mean_<family>(),
# .pair_draw_<family>(), .pair_mean()) take a brmsprep rather than a fit, so
# the tests can call them on a hand-built prep without compiling Stan. With
# keep_observed = FALSE and the same seed, the _y1_given_y2 functions return
# the draws and expectations that posterior_predict() and posterior_epred()
# returned up to 0.10.1, so that results computed with 0.10.1 can be
# reproduced.

.BICOUNT_FAMILIES <- c("bipois", "binegbin")

# --------------------------------------------------------------------------
# Draws and expectations of y1 conditional on y2
# --------------------------------------------------------------------------

.impute_draw_bipois <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y2 <- prep$data$vint1[i]
  # Conditioning a sum of independent Poissons on its total gives a Binomial:
  # N_shared | y2 ~ Binomial(y2, mu / (mu + lambdatwo)). N1 is independent of
  # y2, so it is drawn from its own distribution and added.
  p_shared <- mu / (mu + lambdatwo)
  n_shared <- stats::rbinom(length(mu), size = y2, prob = p_shared)
  n1       <- stats::rpois(length(mu), lambdaone)
  n_shared + n1
}

.impute_mean_bipois <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")   # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  y2 <- prep$data$vint1
  # E[y1 | y2] = y2 * mu / (mu + lambdatwo) + lambdaone, the mean of the
  # Binomial split drawn by .impute_draw_bipois().
  y2_mat <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}

.impute_draw_binegbin <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y2 <- prep$data$vint1[i]
  # A sum of independent negative binomials conditioned on its total is not
  # Binomial. The conditional is
  #   P(N_shared = k | y2) proportional to
  #     NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo),  k = 0..y2,
  # in which shapextwo enters because the count conditioned on is y2. Draw k
  # from P(N_shared = k | y2), then add N1 ~ NB2(lambdaone, shapexone).
  S <- length(mu)
  out <- integer(S)
  for (s in seq_len(S)) {
    k <- 0:y2
    lw <- stats::dnbinom(k,      size = shapes[s],    mu = mu[s],        log = TRUE) +
          stats::dnbinom(y2 - k, size = shapextwo[s], mu = lambdatwo[s], log = TRUE)
    w <- exp(lw - max(lw))
    n_shared <- if (y2 == 0) 0L else sample(k, 1, prob = w)
    out[s] <- n_shared + stats::rnbinom(1, size = shapexone[s], mu = lambdaone[s])
  }
  out
}

.impute_mean_binegbin <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")   # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  shapes    <- brms::get_dpar(prep, "shapes")
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES)
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + lambdaone, with E[N_shared | y2] summed
  # over the same weights that .impute_draw_binegbin() samples from. See
  # .e_shared_given_y2_nb() in utils.R.
  .e_shared_given_y2_nb(mu, lambdatwo, shapes, shapextwo, y2) + lambdaone
}

# --------------------------------------------------------------------------
# Marginal: the pair (y1, y2) drawn afresh
# --------------------------------------------------------------------------

# .draw_components() draws n values of each of the three latent components,
# Poisson when shapes is NULL and NB2 otherwise, in the order N_shared, N1, N2.
# rbicount() and both .pair_draw_<family>() call it, so that simulated data and
# posterior predictive draws of the pair come from one piece of code. The
# caller checks that the three dispersions are all NULL or all supplied.
.draw_components <- function(n, mu, lambdaone, lambdatwo,
                             shapes = NULL, shapexone = NULL,
                             shapextwo = NULL) {
  if (is.null(shapes)) {
    list(
      Ns = stats::rpois(n, mu),
      N1 = stats::rpois(n, lambdaone),
      N2 = stats::rpois(n, lambdatwo)
    )
  } else {
    list(
      Ns = stats::rnbinom(n, size = shapes,    mu = mu),
      N1 = stats::rnbinom(n, size = shapexone, mu = lambdaone),
      N2 = stats::rnbinom(n, size = shapextwo, mu = lambdatwo)
    )
  }
}

# Each .pair_draw_<family>() returns an ndraws x 2 integer matrix for
# observation i, with columns y1 and y2. The shared component is drawn once and
# added to both counts, so Cov(y1, y2) = Var(N_shared).
.pair_draw_bipois <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  comp <- .draw_components(length(mu), mu, lambdaone, lambdatwo)
  cbind(comp$Ns + comp$N1, comp$Ns + comp$N2)
}

.pair_draw_binegbin <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  comp <- .draw_components(
    length(mu),
    mu,
    lambdaone,
    lambdatwo,
    shapes    = shapes,
    shapexone = shapexone,
    shapextwo = shapextwo
  )
  cbind(comp$Ns + comp$N1, comp$Ns + comp$N2)
}

# The marginal means are the same for both families: a negative-binomial
# component with mean m has mean m, whatever its dispersion.
.pair_mean <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  list(y1 = mu + lambdaone, y2 = mu + lambdatwo)
}

# --------------------------------------------------------------------------
# Prep-level drivers
# --------------------------------------------------------------------------

# TRUE on rows where y1 was recorded. A prep from bipois() or binegbin() has no
# vint2, and every row of it is paired. See .y1_obs_at() in utils.R.
.paired_rows <- function(prep) {
  if (is.null(prep$data$vint2)) {
    rep(TRUE, prep$nobs)
  } else {
    prep$data$vint2 == 1L
  }
}

.check_family <- function(family) {
  if (!isTRUE(family %in% .BICOUNT_FAMILIES)) {
    stop(
      "`fit` must use a family from bicountbrms (",
      paste0("`", .BICOUNT_FAMILIES, "`", collapse = " or "),
      "); its family is `", format(family), "`.",
      call. = FALSE
    )
  }
  family
}

.impute_y1_prep <- function(prep, family, type, keep_observed) {
  family <- .check_family(family)
  if (type == "draw") {
    fun <- get(paste0(".impute_draw_", family), mode = "function")
    # Observations are drawn in order, as brms::posterior_predict() draws
    # them, so that the random-number stream matches that of
    # posterior_predict() in 0.10.1.
    out <- do.call(cbind, lapply(seq_len(prep$nobs), fun, prep = prep))
  } else {
    fun <- get(paste0(".impute_mean_", family), mode = "function")
    out <- fun(prep)
  }
  out <- matrix(out, nrow = prep$ndraws, ncol = prep$nobs)
  if (keep_observed) {
    # Draws are made on paired rows too and then replaced, so that the draws
    # on unpaired rows do not depend on keep_observed.
    paired <- which(.paired_rows(prep))
    if (anyNA(prep$data$Y[paired])) {
      stop(
        "`y1` is missing on ", sum(is.na(prep$data$Y[paired])),
        " row(s) indicated as paired. Supply the recorded `y1` in `newdata`, ",
        "or set `keep_observed = FALSE`.",
        call. = FALSE
      )
    }
    if (length(paired)) {
      out[, paired] <- matrix(
        prep$data$Y[paired],
        nrow  = prep$ndraws,
        ncol  = length(paired),
        byrow = TRUE
      )
    }
  }
  out
}

.joint_prep <- function(prep, family, type) {
  family <- .check_family(family)
  S <- prep$ndraws
  N <- prep$nobs
  dn <- list(NULL, NULL, c("y1", "y2"))
  if (type == "draw") {
    fun <- get(paste0(".pair_draw_", family), mode = "function")
    out <- array(NA_integer_, dim = c(S, N, 2L), dimnames = dn)
    for (i in seq_len(N)) {
      out[, i, ] <- fun(i, prep)
    }
  } else {
    m <- .pair_mean(prep)
    out <- array(NA_real_, dim = c(S, N, 2L), dimnames = dn)
    out[, , 1L] <- m$y1
    out[, , 2L] <- m$y2
  }
  out
}

.bicount_prep <- function(object, newdata, re_formula, ndraws, draw_ids, ...) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.", call. = FALSE)
  }
  .check_family(object$family$name)
  # As in brms::posterior_predict(): restructure() updates a fit stored by an
  # older brms, and check_response = FALSE lets `newdata` omit y1, which is
  # then NA in prep$data$Y.
  object <- brms::restructure(object)
  prep <- brms::prepare_predictions(
    object,
    newdata        = newdata,
    re_formula     = re_formula,
    ndraws         = ndraws,
    draw_ids       = draw_ids,
    check_response = FALSE,
    ...
  )
  # Evaluate each dpar once, on the response scale, as
  # brms::posterior_predict() does before its loop over observations;
  # brms::get_dpar(prep, dpar, i) then slices a column rather than
  # re-evaluating the linear predictor for each observation. This loop draws no
  # random numbers.
  for (dp in names(prep$dpars)) {
    prep$dpars[[dp]] <- brms::get_dpar(prep, dpar = dp)
  }
  prep
}

.check_keep_observed <- function(keep_observed) {
  if (!(is.logical(keep_observed) && length(keep_observed) == 1L &&
        !is.na(keep_observed))) {
    stop("`keep_observed` must be TRUE or FALSE.", call. = FALSE)
  }
  keep_observed
}

# --------------------------------------------------------------------------
# Exported functions
# --------------------------------------------------------------------------

#' Draws and expectations of both counts from the joint model
#'
#' @description
#' `posterior_predict_y1_and_y2()` draws the pair `(y1, y2)` from the joint
#' model of the fit, with new shared and source-exclusive components on every
#' draw. `posterior_epred_y1_and_y2()` returns the marginal means of the two
#' counts, `mu + lambdaone` and `mu + lambdatwo`. Neither function conditions
#' on a recorded count. The `y1` slice of `posterior_predict_y1_and_y2()` has
#' the distribution of the draws of `posterior_predict()`, and the `y1` slice
#' of `posterior_epred_y1_and_y2()` equals `posterior_epred()`. Both functions
#' return an array in the shape that brms uses for a multivariate model:
#' posterior draws by observations by response.
#'
#' The two counts are drawn together: the `y1` and `y2` of each draw share one
#' draw of the shared component `N_s`, so the covariance of the drawn counts is
#' `Var(N_s)`, as in the model. The package has no separate function for `y2`,
#' because draws of `y2` made apart from the draws of `y1` would not share a
#' draw of `N_s` with them. Marginal draws of `y2` are `draws[, , "y2"]`, where
#' `draws` is the result of `posterior_predict_y1_and_y2()`.
#'
#' @details
#' For each posterior draw of the parameters, the three components `N_s`, `N1`,
#' and `N2` are drawn independently from their distributions, and
#' `posterior_predict_y1_and_y2()` returns `y1 = N_s + N1` and
#' `y2 = N_s + N2`. The marginal moments are
#'
#' \deqn{E[y_1] = \mu + \lambda_1, \qquad E[y_2] = \mu + \lambda_2,
#'   \qquad \mathrm{Cov}(y_1, y_2) = \mathrm{Var}(N_s),}{E[y1] = mu + lambda1,
#'   E[y2] = mu + lambda2, Cov(y1, y2) = Var(N_s),}
#'
#' with \eqn{\mathrm{Var}(N_s) = \mu}{Var(N_s) = mu} for `bipois` and
#' \eqn{\mathrm{Var}(N_s) = \mu + \mu^2 / \phi_s}{Var(N_s) = mu + mu^2 / phi_s}
#' for `binegbin`, where \eqn{\phi_s}{phi_s} is `shapes`. The variance of each
#' count is the sum of the variances of its two components, for example
#' \eqn{\mathrm{Var}(y_1) = \mu + \mu^2/\phi_s + \lambda_1 +
#' \lambda_1^2/\phi_{x1}}{Var(y1) = mu + mu^2/phi_s + lambda1 +
#' lambda1^2/phi_x1} for `binegbin`.
#'
#' A predictive check of the pair, such as of the variance of the difference
#' `y1 - y2`, needs draws of both counts from `posterior_predict_y1_and_y2()`.
#' Subtracting the recorded `y2` from draws of `posterior_predict()` treats the
#' two counts as independent and overstates that variance by
#' \eqn{2\,\mathrm{Var}(N_s)}{2 Var(N_s)}.
#'
#' **New data.** The recorded `y1`, the recorded `y2`, and the observation indicator
#' do not affect the result. The column `y1` may be omitted from `newdata`, but
#' brms requires every variable named in the `vint()` term, so `newdata` must
#' contain `y2` and, for a fit from a `_partialobs()` constructor, the
#' observation indicator.
#' The values in these columns are not used, so any integers will do, such as
#' `0` or `-1`. The same is true for `newdata` in `posterior_predict()` and
#' `posterior_epred()`.
#'
#' @inheritParams posterior_predict_y1_given_y2
#'
#' @return A numeric array with dimensions posterior draws x observations x 2,
#'   with the third dimension named `y1` and `y2`. The values from
#'   `posterior_predict_y1_and_y2()` are integers.
#'
#' @seealso [posterior_predict_y1_given_y2()] for draws of `y1` conditional
#'   on the recorded `y2`; `vignette("prediction")` for a summary of the
#'   prediction functions.
#'
#' @examples
#' \dontrun{
#' library(brms)
#'
#' # 500 units. y2 is recorded on every unit, and y1 on about 20% of units.
#' set.seed(20261115)
#' dat <- rbicount(
#'   500,
#'   mu        = 4,
#'   lambdaone = 1,
#'   lambdatwo = 2,
#'   shapes    = 2,
#'   shapexone = 4,
#'   shapextwo = 4
#' )
#' dat$y1_obs  <- rbinom(nrow(dat), 1, 0.2)
#' dat$y1_true <- dat$y1
#' dat$y1[dat$y1_obs == 0] <- -1L
#'
#' fit <- brm(
#'   bf_mfd(
#'     y1 | vint(y2, y1_obs),
#'     logM      ~ 1,
#'     logitf    ~ 1,
#'     delta     ~ 1,
#'     kappas    ~ 1,
#'     kappaxone ~ 1,
#'     kappaxtwo ~ 1
#'   ),
#'   family   = binegbin_partialobs(),
#'   stanvars = binegbin_partialobs_stanvars(),
#'   data     = dat,
#'   prior    = c(
#'     prior(normal(1.5, 1), class = "b", nlpar = "logM"),
#'     prior(normal(0, 1.5), class = "b", nlpar = "logitf"),
#'     prior(normal(0, 0.5), class = "b", nlpar = "delta"),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappas",    lb = 0),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappaxone", lb = 0),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappaxtwo", lb = 0)
#'   )
#' )
#'
#' # Predictive check of the variance of the difference, on the paired rows
#' paired <- dat[dat$y1_obs == 1, ]
#' yrep <- posterior_predict_y1_and_y2(fit, newdata = paired, ndraws = 500)
#' var_d_rep <- apply(yrep[, , "y1"] - yrep[, , "y2"], 1, var)
#' hist(var_d_rep)
#' abline(v = var(paired$y1 - paired$y2))
#'
#' # Expected counts of the two sources on each paired unit
#' ep <- posterior_epred_y1_and_y2(fit, newdata = paired)
#' dim(ep)
#' }
#' @export
posterior_predict_y1_and_y2 <- function(object,
                                        newdata    = NULL,
                                        re_formula = NULL,
                                        ndraws     = NULL,
                                        draw_ids   = NULL,
                                        ...) {
  prep <- .bicount_prep(
    object,
    newdata    = newdata,
    re_formula = re_formula,
    ndraws     = ndraws,
    draw_ids   = draw_ids,
    ...
  )
  .joint_prep(
    prep,
    family = prep$family$name,
    type   = "draw"
  )
}

#' @rdname posterior_predict_y1_and_y2
#' @export
posterior_epred_y1_and_y2 <- function(object,
                                      newdata    = NULL,
                                      re_formula = NULL,
                                      ndraws     = NULL,
                                      draw_ids   = NULL,
                                      ...) {
  prep <- .bicount_prep(
    object,
    newdata    = newdata,
    re_formula = re_formula,
    ndraws     = ndraws,
    draw_ids   = draw_ids,
    ...
  )
  .joint_prep(
    prep,
    family = prep$family$name,
    type   = "mean"
  )
}

#' Draws and expectations of y1 given y2
#'
#' @description
#' `posterior_predict_y1_given_y2()` draws `y1` from its distribution
#' conditional on the recorded `y2` and the model
#' parameters, \eqn{p(y_1 \mid y_2, \theta)}{p(y1 | y2, theta)}.
#' `posterior_epred_y1_given_y2()` returns the conditional expectation
#' \eqn{E[y_1 \mid y_2]}{E[y1 | y2]}. The main use of both functions is to
#' impute `y1` on the rows where `y1` was not recorded. On an
#' unpaired row of a fit from [bipois_partialobs()] or
#' [binegbin_partialobs()], the draws are the imputed values of the
#' unrecorded `y1`, and the expectation is the mean of those imputed
#' values.
#'
#' `posterior_predict()` and `posterior_epred()` do not condition on `y2`.
#' For a fit from any of the four constructors, these two functions return
#' draws from the marginal distribution of `y1` and its marginal mean
#' `mu + lambdaone`, as for any other brms family.
#' [posterior_predict_y1_and_y2()] draws both counts.
#'
#' @details
#' Both families construct the pair as `y1 = N_s + N1` and `y2 = N_s + N2`,
#' where the shared component `N_s` has mean `mu` and dispersion `shapes`, and
#' the source-exclusive components `N1` and `N2` have means `lambdaone` and
#' `lambdatwo` and dispersions `shapexone` and `shapextwo`. The three
#' components are independent given the parameters, so `N1` is independent of
#' `y2`, and \eqn{y_1 \mid y_2}{y1 | y2} is the sum of
#' \eqn{N_s \mid y_2}{N_s | y2} and an independent draw of `N1`. The
#' conditional distribution of the shared component is
#'
#' \deqn{P(N_s = k \mid y_2) \propto p(k \mid \mu, \phi_s)\,
#'   p(y_2 - k \mid \lambda_2, \phi_{x2}), \quad k = 0, \ldots, y_2,}{P(N_s =
#'   k | y2) proportional to p(k | mu, phi_s) p(y2 - k | lambda2, phi_x2),
#'   k = 0, ..., y2,}
#'
#' where \eqn{p(\cdot \mid m, \phi)}{p(. | m, phi)} is the negative-binomial
#' probability mass function with mean \eqn{m}{m} and inverse dispersion
#' \eqn{\phi}{phi} for `binegbin`, and the Poisson probability mass function
#' with mean \eqn{m}{m} for `bipois`. For `bipois`, this distribution is
#' Binomial with size \eqn{y_2}{y2} and probability
#' \eqn{\mu / (\mu + \lambda_2)}{mu / (mu + lambda2)}. For `binegbin`, the
#' distribution has no closed form and is evaluated over the
#' \eqn{y_2 + 1}{y2 + 1} values of \eqn{k}{k}. In both families,
#'
#' \deqn{E[y_1 \mid y_2] = E[N_s \mid y_2] + \lambda_1.}{E[y1 | y2] =
#'   E[N_s | y2] + lambda1.}
#'
#' The exclusive-component dispersion in the conditional is `shapextwo`, the
#' dispersion of the second exclusive component, because the count conditioned
#' on is `y2`.
#'
#' **Recorded values of `y1`.** With `keep_observed = FALSE`, the default, every
#' row receives a conditional draw or the conditional expectation, including
#' the paired rows. On a paired row, the conditional draw is the quantity to
#' compare with the recorded `y1` in a conditional predictive check. With
#' `keep_observed = TRUE`, every draw on a paired row is the recorded `y1`, and
#' the draws on unpaired rows are unchanged. Summed over rows, the draws
#' returned with `keep_observed = TRUE` are draws of the total of `y1`:
#' recorded where `y1` was recorded, and imputed where `y1` was not recorded.
#' To reproduce results from `posterior_predict()` and `posterior_epred()` in
#' bicountbrms 0.10.1, use these two functions with `keep_observed = FALSE`;
#' given the same seed, the draws are identical.
#'
#' **New data.** `newdata` must contain the recorded `y2` of each row, because
#' the result depends on `y2`. For a fit from a `_partialobs()` constructor,
#' brms also requires the observation indicator named in `vint()`. The column `y1`
#' may be omitted, except that with `keep_observed = TRUE` a row indicated as
#' paired needs its recorded `y1`, and the function stops with an error if
#' that `y1` is missing.
#'
#' Stored fits whose family object declares the pre-0.7.0 rate names
#' (`lambdaem`, `lambdalb`) or a single exclusive-component dispersion
#' (`shapex`) are read under those names, as by the brms methods.
#'
#' @param object A `brmsfit` from [bipois()], [bipois_partialobs()],
#'   [binegbin()], or [binegbin_partialobs()].
#' @param newdata An optional data frame; if `NULL`, the data of the fit are
#'   used. See Details.
#' @param re_formula A formula for the group-level effects to include, passed
#'   to [brms::prepare_predictions()]. `NULL` (the default) includes all
#'   group-level effects, and `NA` includes none.
#' @param ndraws,draw_ids Passed to [brms::prepare_predictions()] to select
#'   posterior draws.
#' @param keep_observed Logical. If `FALSE` (the default), every row returns
#'   a conditional draw or the conditional expectation. If `TRUE`, paired rows
#'   return the recorded `y1`.
#' @param ... Further arguments to [brms::prepare_predictions()], such as
#'   `allow_new_levels`.
#'
#' @return A numeric matrix with posterior draws in rows and observations in
#'   columns. The values from `posterior_predict_y1_given_y2()` are integers.
#'
#' @seealso [posterior_predict_y1_and_y2()] for draws of both counts from the
#'   marginal model; `vignette("prediction")` for a summary of the prediction
#'   functions.
#'
#' @examples
#' \dontrun{
#' library(brms)
#'
#' # 500 units. y2 is recorded on every unit, and y1 on about 20% of units.
#' set.seed(20261115)
#' dat <- rbicount(
#'   500,
#'   mu        = 4,
#'   lambdaone = 1,
#'   lambdatwo = 2,
#'   shapes    = 2,
#'   shapexone = 4,
#'   shapextwo = 4
#' )
#' dat$y1_obs  <- rbinom(nrow(dat), 1, 0.2)
#' dat$y1_true <- dat$y1
#' dat$y1[dat$y1_obs == 0] <- -1L
#'
#' fit <- brm(
#'   bf_mfd(
#'     y1 | vint(y2, y1_obs),
#'     logM      ~ 1,
#'     logitf    ~ 1,
#'     delta     ~ 1,
#'     kappas    ~ 1,
#'     kappaxone ~ 1,
#'     kappaxtwo ~ 1
#'   ),
#'   family   = binegbin_partialobs(),
#'   stanvars = binegbin_partialobs_stanvars(),
#'   data     = dat,
#'   prior    = c(
#'     prior(normal(1.5, 1), class = "b", nlpar = "logM"),
#'     prior(normal(0, 1.5), class = "b", nlpar = "logitf"),
#'     prior(normal(0, 0.5), class = "b", nlpar = "delta"),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappas",    lb = 0),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappaxone", lb = 0),
#'     prior(normal(0, 1),   class = "b", nlpar = "kappaxtwo", lb = 0)
#'   )
#' )
#'
#' # Draws of y1 on every unit: recorded where y1_obs is 1,
#' # imputed from the y2 of each unit where y1_obs is 0.
#' y1_draws <- posterior_predict_y1_given_y2(fit, keep_observed = TRUE)
#'
#' # Posterior of the source-1 total, against the simulated total
#' quantile(rowSums(y1_draws), probs = c(0.025, 0.5, 0.975))
#' sum(dat$y1_true)
#'
#' # Imputed draws and E[y1 | y2] on the units where source 1 was not
#' # recorded only
#' unrec <- dat$y1_obs == 0
#' y1_imp <- posterior_predict_y1_given_y2(fit)[, unrec]
#' ep_imp <- posterior_epred_y1_given_y2(fit)[, unrec]
#' }
#' @export
posterior_predict_y1_given_y2 <- function(object,
                                          newdata       = NULL,
                                          re_formula    = NULL,
                                          ndraws        = NULL,
                                          draw_ids      = NULL,
                                          keep_observed = FALSE,
                                          ...) {
  .check_keep_observed(keep_observed)
  prep <- .bicount_prep(
    object,
    newdata    = newdata,
    re_formula = re_formula,
    ndraws     = ndraws,
    draw_ids   = draw_ids,
    ...
  )
  .impute_y1_prep(
    prep,
    family        = prep$family$name,
    type          = "draw",
    keep_observed = keep_observed
  )
}

#' @rdname posterior_predict_y1_given_y2
#' @export
posterior_epred_y1_given_y2 <- function(object,
                                        newdata       = NULL,
                                        re_formula    = NULL,
                                        ndraws        = NULL,
                                        draw_ids      = NULL,
                                        keep_observed = FALSE,
                                        ...) {
  .check_keep_observed(keep_observed)
  prep <- .bicount_prep(
    object,
    newdata    = newdata,
    re_formula = re_formula,
    ndraws     = ndraws,
    draw_ids   = draw_ids,
    ...
  )
  .impute_y1_prep(
    prep,
    family        = prep$family$name,
    type          = "mean",
    keep_observed = keep_observed
  )
}
