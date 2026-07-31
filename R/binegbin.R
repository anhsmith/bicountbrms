# ==========================================================================
# binegbin: joint bivariate Negative-Binomial via trivariate reduction
#
# The overdispersed sibling of bipois (see bipois.R). Same trivariate-
# reduction construction -- y1 = N_shared + N1, y2 = N_shared + N2,
# the three latent counts mutually independent given their rates -- but each
# latent count is Negative-Binomial rather than Poisson:
#
#   N_shared ~ NB2(mu,        shapes)   (shared component; drives correlation)
#   N1       ~ NB2(lambdaone, shapex)   (source-1-only excess)
#   N2       ~ NB2(lambdatwo, shapex)   (source-2-only excess)
#
# NB2(m, phi) is Stan's neg_binomial_2 / R's dnbinom(size = phi, mu = m):
# mean m, variance m + m^2/phi. shapes is the shared-component dispersion and
# shapex the excess dispersion, shared by N1 and N2 -- two dispersion
# parameters rather than three, since the two private components are not
# separately identified in the data these families are built for.
#
# WHY NEGBIN AND NOT AN OLRE ON bipois. The plain-Poisson bipois cannot be
# overdispersed (Var == mean for each latent count), so it underfits the
# real marginal variances badly -- in one motivating dataset by ~10x
# (Var(y1) fitted 16.6 against 179 observed) and Var(d) by ~3.5x. The
# obvious fix -- add a per-set
# observation-level random effect (OLRE) on the excess components -- FAILS
# synthetic recovery: with one bivariate observation per set but three
# per-set latent deviates (mu-OLRE + two excess OLREs), the excess deviates
# act as residual-absorbers, their population SD collapses toward the prior
# mode, and drawing fresh deviates does NOT regenerate the observed spread
# (recovered excess SD 0.37 vs true 0.85; fresh-deviate Var(d) 2.9 vs true
# 19.2). A conditional posterior-predictive check hides this completely --
# only a marginal (fresh-deviate) check exposes it. binegbin carries the
# dispersion in SCALAR shapes/shapex instead, estimated from aggregate
# mean-variance mismatch across sets -- identifiable, no per-set overfitting,
# clean marginal PPC, and consistent with the review-track NegBin models.
#
# LIKELIHOOD. N_shared is unobserved and marginalised out analytically,
# exactly as in bipois -- the sum structure is identical, only the component
# pmfs change from Poisson to NegBin:
#
#   P(y1=x, y2=y) = sum_{k=0}^{min(x,y)}
#     NB2(k | mu, shapes) NB2(x-k | lambdaone, shapex) NB2(y-k | lambdatwo, shapex)
#
# This is NOT a "two stacked marginalisations" problem
# (Gamma-mixing a Poisson while keeping it Poisson): the components are
# DIRECTLY NegBin, so the marginalisation sum is the same finite sum as
# bipois with neg_binomial_2_lpmf swapped in for poisson_lpmf. The bipois
# incremental recurrence does not carry over cleanly (NegBin consecutive-term
# ratios are less simple than Poisson's), so binegbin_lpmf evaluates the sum
# directly via log_sum_exp. m = min(y1, y2) is bounded by the data
# (~50-60 terms at most for this project), so the direct sum is not a
# performance concern -- the same reasoning bipois's own docs give for why no
# large-argument branch is needed.
#
# Validation -- grid cross-check of the Stan lpmf against the independent R
# brute-force reference to ~1e-14, normalisation to 1, the moment identities,
# the Poisson-limit reduction to bipois, and end-to-end parameter recovery
# with a coverage assessment -- is in tests/testthat/test-binegbin.R.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Joint bivariate-Negative-Binomial custom family for brms
#'
#' @description
#' Overdispersed sibling of [bipois()]. Returns a brms custom family for the
#' joint distribution of a matched count pair `(y1, y2)` via trivariate
#' reduction with Negative-Binomial (rather than Poisson) latent components:
#' `y1 = N_shared + N1`, `y2 = N_shared + N2`, with
#' `N_shared ~ NB2(mu, shapes)`, `N1 ~ NB2(lambdaone, shapex)`,
#' `N2 ~ NB2(lambdatwo, shapex)` mutually independent given their rates.
#' `NB2(m, phi)` has mean `m` and variance `m + m^2/phi` (Stan
#' `neg_binomial_2`; R `dnbinom(size = phi, mu = m)`).
#'
#' Five dpars: the three rates (`mu` = shared rate, `lambdaone`/`lambdatwo` =
#' the two source-specific rates) plus two dispersions -- `shapes` for the
#' shared component and `shapex` shared across the two excess components. All
#' five use `link = "log"`. Supply the excess rates through a non-linear
#' formula without an explicit `exp()` (the log link applies it):
#' `nlf(lambdaone ~ lamx)` gives `lambdaone = exp(lamx)`.
#'
#' See the `binegbin.R` file header for why Negative-Binomial components are
#' used instead of an observation-level random effect on [bipois()] -- briefly,
#' the random-effect version fails synthetic recovery, because with one
#' observed pair but three latent deviates per unit the excess deviates act as
#' residual absorbers.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2) ~ 1,
#'        mu ~ 1 + (1 | vessel),
#'        nlf(lambdaone ~ lamx), nlf(lambdatwo ~ lamx), lamx ~ 1,
#'        shapes ~ 1, shapex ~ 1, nl = TRUE),
#'     family   = binegbin(),
#'     stanvars = binegbin_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Forced `mu` naming, and `y2` via `vint()`.** Identical conventions to
#' [bipois()] -- `mu` is brms's mandatory dpar name, here bound to the shared
#' component's rate (`lambda_shared`), not a mean of either response; `y2`
#' travels as supplementary integer data through `vint()` because
#' `custom_family()` declares a single response column. See [bipois()] for
#' the full explanation, including why the rates are spelled `lambdaone`/
#' `lambdatwo` in code but written \eqn{\lambda_1}{lambda_1}/
#' \eqn{\lambda_2}{lambda_2} in the documentation.
#'
#' **Order of dpars matters for the generated Stan call.** brms generates
#' `target += binegbin_lpmf(Y[n] | mu[n], lambdaone[n], lambdatwo[n],
#' shapes[n], shapex[n], vint1[n])` -- dpars in the order declared here, then
#' vint args. `binegbin_stan_funs` declares `binegbin_lpmf` with exactly this
#' signature; reordering one without the other silently swaps which rate or
#' dispersion governs which component.
#'
#' @return A brms custom_family object.
#' @export
binegbin <- function() {
  brms::custom_family(
    name  = "binegbin",
    dpars = c("mu", "lambdaone", "lambdatwo", "shapes", "shapex"),
    links = c("log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0),
    type  = "int",
    vars  = "vint1[n]"
  )
}

#' @rdname binegbin
#' @export
binegbin_stanvars <- function() {
  brms::stanvar(block = "functions", scode = binegbin_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Direct marginalisation sum (see file header for why not a recurrence).
# Every term is a sum of three neg_binomial_2 log-densities; accumulated via
# log_sum_exp over k = 0..min(y1, y2). neg_binomial_2_lpmf(0 | m, phi) is
# well-defined, so the k = 0 term and zero-count responses need no special
# casing.
binegbin_stan_funs <- "
  real binegbin_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
                     real shapes, real shapex, int y2) {
    int m = min(y1, y2);
    vector[m + 1] lp;
    for (k in 0:m) {
      lp[k + 1] = neg_binomial_2_lpmf(k      | mu,        shapes)
                + neg_binomial_2_lpmf(y1 - k | lambdaone, shapex)
                + neg_binomial_2_lpmf(y2 - k | lambdatwo, shapex);
    }
    return log_sum_exp(lp);
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Independent brute-force evaluation of the same sum via R's dnbinom (an
# independent route from Stan's neg_binomial_2), used to validate the Stan
# lpmf and to power log_lik_binegbin()/posterior_epred_binegbin() post-hoc.
# Internal reference only, mirroring bipois_lpmf_r's role.
binegbin_lpmf_r <- function(y1, y2, mu, lambdaone, lambdatwo, shapes, shapex) {
  n <- max(length(y1), length(y2), length(mu), length(lambdaone),
           length(lambdatwo), length(shapes), length(shapex))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)
  shapes    <- rep_len(shapes, n)
  shapex    <- rep_len(shapex, n)

  vapply(seq_len(n), function(i) {
    m <- min(y1[i], y2[i])
    k <- 0:m
    log_terms <- stats::dnbinom(k,        size = shapes[i], mu = mu[i],        log = TRUE) +
      stats::dnbinom(y1[i] - k,           size = shapex[i], mu = lambdaone[i], log = TRUE) +
      stats::dnbinom(y2[i] - k,           size = shapex[i], mu = lambdatwo[i], log = TRUE)
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname binegbin
#' @export
#' @keywords internal
log_lik_binegbin <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapex    <- brms::get_dpar(prep, "shapex", i = i)
  y1 <- prep$data$Y[i]
  y2 <- prep$data$vint1[i]
  binegbin_lpmf_r(y1, y2, mu, lambdaone, lambdatwo, shapes, shapex)
}

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_predict_binegbin <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapex    <- brms::get_dpar(prep, "shapex", i = i)
  y2 <- prep$data$vint1[i]
  # y1 predicted conditional on the real observed y2. Unlike bipois, the
  # conditional split N_shared | y2 is NOT Binomial (a NegBin sum condition
  # is not Binomial); it is P(N_shared = k | y2) proportional to
  # NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex) over k = 0..y2.
  # Sample that discrete conditional, then add a fresh N1 ~ NB2(lambdaone,
  # shapex).
  S <- length(mu)
  out <- integer(S)
  for (s in seq_len(S)) {
    k <- 0:y2
    lw <- stats::dnbinom(k,      size = shapes[s], mu = mu[s],        log = TRUE) +
          stats::dnbinom(y2 - k, size = shapex[s], mu = lambdatwo[s], log = TRUE)
    w <- exp(lw - max(lw))
    n_shared <- if (y2 == 0) 0L else sample(k, 1, prob = w)
    out[s] <- n_shared + stats::rnbinom(1, size = shapex[s], mu = lambdaone[s])
  }
  out
}

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_epred_binegbin <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + lambdaone. For binegbin there is no
  # closed-form E[N_shared | y2] as clean as bipois's y2 * mu/(mu+lambdatwo);
  # this returns the analogous point approximation using the marginal shared
  # fraction mu/(mu+lambdatwo), adequate for epred display (posterior_predict
  # uses the exact discrete conditional).
  y2_mat <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}
