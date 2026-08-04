# ==========================================================================
# bipois: joint bivariate Poisson via trivariate reduction
# (Holgate 1964; recurrence: "A fast way to calculate the bivariate poisson
# in STAN", stan-users Google Group, March 2016 -- thread SUcp-ktkXn4,
# posted by Andre, refined by Bob Carpenter -- adapted below from
# log-rate/`poisson_log_lpmf` parameterisation to the natural-scale-rate
# parameterisation this package's dpars use post-link.)
#
# Unlike a difference family -- the Skellam, discrete Laplace and discrete
# normal families in the companion package skellambrms
# (https://github.com/anhsmith/skellambrms), all of which model a single
# response d = y1 - y2 -- bipois models the *joint* pair (y1, y2)
# directly. That avoids regressing the difference on one of its own
# components (a `d ~ y2` design), which induces regression to the mean: y2
# appears on both sides, so the fitted slope is biased towards -1 by the
# shared sampling error alone. Both y1 and y2 are generated as
# y1 = N_shared + N1, y2 = N_shared + N2, with N_shared ~ Poisson(mu),
# N1 ~ Poisson(lambdaone), N2 ~ Poisson(lambdatwo) mutually independent given
# their rates. N_shared is unobserved and marginalised out analytically:
#
#   P(y1=x, y2=y) = sum_{k=0}^{min(x,y)}
#     P(N_shared=k) P(N1=x-k) P(N2=y-k)
#
# rather than evaluated as a naive per-k sum of three separate Poisson
# lpmfs, the Stan implementation below uses the cited bipois2 incremental
# recurrence: the k=0 term is computed directly, and each subsequent term
# is obtained from the previous one by a ratio update in log space,
# accumulated via log_sum_exp -- see bipois_stan_funs below for the
# derivation of that ratio and its correspondence to the cited source.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Joint bivariate-Poisson custom family for brms
#'
#' @description
#' Returns a brms custom family for the joint distribution of a matched pair
#' of counts, `(y1, y2)`, constructed via trivariate reduction: `y1 =
#' N_shared + N1`, `y2 = N_shared + N2`, with `N_shared ~
#' Poisson(mu)`, `N1 ~ Poisson(lambdaone)`, `N2 ~ Poisson(lambdatwo)`
#' mutually independent given their rates. All three rates are link = "log".
#' Modelling the pair jointly avoids regressing the difference on one of its
#' own components (a `d = y1 - y2 ~ y2` design), which induces regression to
#' the mean.
#'
#' `y1` is the family's response; `y2` is passed in as supplementary
#' integer data via brms's `vint()` addition term, since brms's
#' `custom_family()` machinery is built around a single declared response
#' column -- see Details.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2) ~ ...),
#'     family   = bipois(),
#'     stanvars = bipois_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Naming note.** `brms::custom_family()` requires a dpar literally named
#' `"mu"` (`stop2("All families must have a 'mu' parameter.")`,
#' unconditional), whatever the family actually calls that quantity. Here it is
#' bound to `lambda_shared`, the rate of the
#' component shared between `y1` and `y2` -- not a mean of either
#' response individually. `lambdaone` (source-1-only rate) and `lambdatwo`
#' (source-2-only rate) are the other two dpars, plainly named (no forced
#' reinterpretation needed for those two).
#'
#' **Why `lambdaone`, not `lambda1`.** `custom_family()` rejects dpar names
#' ending in a digit (`stop2("'dpars' should not end with a number.")`), as
#' well as dots and underscores. The documentation therefore writes these
#' rates as \eqn{\lambda_1}{lambda_1} and \eqn{\lambda_2}{lambda_2} while the
#' code must spell them `lambdaone`/`lambdatwo`. See the notation table in
#' the package README.
#'
#' **Why `y2` travels via `vint()`, not as a second response.** brms's
#' `custom_family()` API supports exactly one declared response column
#' (`Y`) plus optional supplementary integer/real data (`vint()`/`vreal()`
#' addition terms) -- the same mechanism used for, e.g., binomial trial
#' counts in the brms custom-families vignette. There is no
#' *undeclared-response* concept for a genuinely joint two-count
#' likelihood; `vint(y2)` is the correct fit for that gap, not a
#' workaround. This does mean `y2` is *not* itself treated as
#' brms-modelled response data (no missing-value handling, no
#' resp_*() addition terms apply to it) -- it is fixed, observed
#' per-row data, consistent with the fact that every row used here comes
#' from the matched (both-observed) subset.
#'
#' **Order of dpars matters for the generated Stan call.** brms generates
#' `target += bipois_lpmf(Y[n] | mu[n], lambdaone[n], lambdatwo[n],
#' vint1[n])` -- dpars in the order declared here, then vint/vreal args in
#' the order declared in `vars`. `bipois_stan_funs` (stanfunctions via
#' `bipois_stanvars()`) declares `bipois_lpmf` with exactly this argument
#' order; changing the order here without changing the Stan signature (or
#' vice versa) silently swaps which rate governs which count.
#'
#' @return A brms custom_family object.
#' @export
bipois <- function() {
  brms::custom_family(
    name  = "bipois",
    dpars = c("mu", "lambdaone", "lambdatwo"),  # mu = lambda_shared -- see Details
    links = c("log", "log", "log"),
    lb    = c(0, 0, 0),
    type  = "int",
    vars  = "vint1[n]"
  )
}

#' @rdname bipois
#' @export
bipois_stanvars <- function() {
  brms::stanvar(block = "functions", scode = bipois_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Derivation of the recurrence, in this package's natural-scale-rate terms
# (the cited stan-users thread works in log-rates via poisson_log_lpmf;
# translated here since this package's dpars are already inv-link-
# transformed to natural scale by the time they reach the lpmf, the same
# convention as every other family in this package: binegbin_lpmf likewise
# takes `mu`, not `log(mu)`).
#
# Write lambda_shared = mu, and let term(k) = P(N_shared=k) P(N1=r-k)
# P(N2=s-k) for r = y1, s = y2, m = min(r,s). Then:
#
#   term(0) = exp(-mu) * [lambdaone^r exp(-lambdaone) / r!]
#                       * [lambdatwo^s exp(-lambdatwo) / s!]
#
# so log(term(0)) = poisson_lpmf(r|lambdaone) + poisson_lpmf(s|lambdatwo) - mu
# (the `+lambdaone +lambdatwo` inside each poisson_lpmf cancels algebraically
# against the `-lambdaone -lambdatwo` that would otherwise appear from
# factoring exp(-mu-lambdaone-lambdatwo) out front -- this is exactly the
# same cancellation used in the cited thread's `ss <- poisson_log_log(r,
# mu1) + poisson_log_log(s, mu2) - exp(mu3)` starting term).
#
# The ratio between consecutive terms is:
#   term(k) / term(k-1) = [(r-k+1)/lambdaone] * [(s-k+1)/lambdatwo] * [mu/k]
# (one fewer factor of lambdaone in the N1 count as k increases by one,
# one fewer of lambdatwo, one more of mu -- and the corresponding
# factorial/combinatorial adjustment (r-k+1), (s-k+1), 1/k). In log space:
#   log(term(k)) = log(term(k-1)) + log(r-k+1) + log(s-k+1) - log(k)
#                  + log(mu) - log(lambdaone) - log(lambdatwo)
# which is exactly the cited thread's per-step update (`log_s <- log_s +
# log(r-k+1) + mus + log(s-k+1) - log(k)`, with `mus = -mu1-mu2+mu3`
# translating directly to `log(mu) - log(lambdaone) - log(lambdatwo)` once
# mu1, mu2, mu3 are read as log-rates). Each term is accumulated into the
# running total via log_sum_exp -- no separate normalising-constant term is
# needed, since this sum is finite (m+1 terms) and exactly equals the
# marginalised joint log-likelihood, not a survival function.
#
# No large-argument blowup risk: m = min(y1, y2) is bounded by the data itself
# (a few tens of terms at most for counts of the size these families are built
# for), and every term here is a sum/difference of ordinary Poisson
# log-densities, evaluated at arguments the data bounds. No iteration cap or
# normal-approximation branch is needed for that reason -- in contrast to the
# Bessel-function evaluation a Skellam log-CCDF requires, where the argument is
# unbounded during warmup and a fallback branch is unavoidable (see
# skellambrms).
bipois_stan_funs <- "
  real bipois_lpmf(int y1, real mu, real lambdaone, real lambdatwo, int y2) {
    int m = min(y1, y2);
    real ss = poisson_lpmf(y1 | lambdaone) + poisson_lpmf(y2 | lambdatwo) - mu;
    if (m > 0) {
      real log_ratio = log(mu) - log(lambdaone) - log(lambdatwo);
      real log_term = ss;
      for (k in 1:m) {
        log_term += log(y1 - k + 1) + log(y2 - k + 1) - log(k) + log_ratio;
        ss = log_sum_exp(ss, log_term);
      }
    }
    return ss;
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Direct brute-force evaluation of the same marginalised sum, computed
# term-by-term from the original P(N_shared=k)*P(N1=x-k)*P(N2=y-k)
# definition rather than via the recurrence --
# an independent route to the same quantity, used to validate the Stan
# recurrence above and to power log_lik_bipois()/posterior_epred_bipois()
# (evaluated post-hoc, not inside the sampler's hot loop, so there is no
# reason to use the recurrence's algebraic shortcuts here). Not exported:
# internal reference only, exactly the role binegbin_lpmf_r plays for its own
# family.
bipois_lpmf_r <- function(y1, y2, mu, lambdaone, lambdatwo) {
  n <- max(length(y1), length(y2), length(mu), length(lambdaone), length(lambdatwo))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  vapply(seq_len(n), function(i) {
    m <- min(y1[i], y2[i])
    k <- 0:m
    log_terms <- stats::dpois(k, mu[i], log = TRUE) +
      stats::dpois(y1[i] - k, lambdaone[i], log = TRUE) +
      stats::dpois(y2[i] - k, lambdatwo[i], log = TRUE)
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname bipois
#' @export
#' @keywords internal
log_lik_bipois <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared -- see Details
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y1 <- prep$data$Y[i]
  y2 <- prep$data$vint1[i]
  bipois_lpmf_r(y1, y2, mu, lambdaone, lambdatwo)
}

#' @rdname bipois
#' @export
#' @keywords internal
posterior_predict_bipois <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared -- see Details
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y2 <- prep$data$vint1[i]
  # y2 is fixed, observed data (not itself re-simulated -- see "Why y2
  # travels via vint()" in ?bipois). Consistent with that, y1 is
  # predicted *conditional on the real y2*, via the closed-form conditional
  # split:
  # N_shared | y2 ~ Binomial(y2, mu / (mu + lambdatwo)); N1 fresh from
  # its own marginal; y1 = N_shared + N1.
  p_shared <- mu / (mu + lambdatwo)
  n_shared <- stats::rbinom(length(mu), size = y2, prob = p_shared)
  n1       <- stats::rpois(length(mu), lambdaone)
  n_shared + n1
}

#' @rdname bipois
#' @export
#' @keywords internal
posterior_epred_bipois <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")   # lambda_shared -- see Details
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + E[N1] = y2 * mu/(mu+lambdatwo) + lambdaone,
  # the same conditional split as posterior_predict_bipois above, in
  # expectation rather than simulated.
  y2_mat <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}
