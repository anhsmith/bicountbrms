# ==========================================================================
# bipois_cens: censoring-aware bivariate Poisson
#
# The censoring-aware extension of bipois (see bipois.R), standing in the same
# relation to it as binegbin_cens does to binegbin. Same generative model --
# y1 = N_shared + N1, y2 = N_shared + N2 with
#
#   N_shared ~ Poisson(mu)          (shared component; drives correlation)
#   N1       ~ Poisson(lambdaone)   (source-1-only excess)
#   N2       ~ Poisson(lambdatwo)   (source-2-only excess)
#
# and the same three log-linked dpars -- but each row now carries a second
# supplementary integer, an observation flag y1_obs, alongside y2:
#
#   y1_obs == 1 (matched set):  full joint bipois lpmf on (y1, y2).
#   y1_obs == 0 (y2-only set):  the y2 MARGINAL of the SAME bivariate model.
#
# THE CENSORED BRANCH IS ANALYTIC, AND THAT IS THE POINT. For binegbin_cens
# the y1-integrated marginal is a convolution that must be evaluated as a sum,
#
#   P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo),
#
# because NB2 + NB2 is not NB2 unless the two components share their success
# probability. The Poisson case has no such obstruction: a sum of independent
# Poissons is Poisson, so the same convolution collapses exactly to
#
#   y2 ~ Poisson(mu + lambdatwo),
#
# one closed-form line with no sum, no cutoff and no accumulated rounding. Two
# consequences beyond the family being cheap to evaluate:
#
#   * It supplies an ANALYTIC REFERENCE for binegbin_cens's censored branch,
#     which the package otherwise lacks. That branch's marginal-identity test
#     checks the NB sum against itself -- summing the matched branch over y1
#     and comparing to the censored branch -- so both sides share any error in
#     the convolution. Taking the NB family to its Poisson limit and comparing
#     against this closed form is an independent check of the same code path.
#     See test-bipois_cens.R.
#   * It avoids a boundary. A user whose counts really are equidispersed must
#     otherwise fit binegbin_cens with the dispersions pressed against their
#     Poisson limit (shape -> Inf, i.e. kappa -> 0), which is exactly where
#     sampling degrades. A dedicated family sidesteps that.
#
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS. Unchanged from
# binegbin_cens, and the argument does not depend on the component
# distribution: the unmatched rows never observe y1, so a plain bipois fit
# could only use the matched rows, yet those rows' y2 is a draw from the same
# bivariate model and still informs mu, lambdatwo and any group-level
# structure. Integrating the unobserved margin out pools every row under one
# coherent likelihood. See binegbin_cens.R for the full statement.
#
# WHAT EACH RATE IS IDENTIFIED FROM. lambdatwo enters both branches, so every
# row sharpens it. lambdaone enters only the matched branch and is identified
# solely by the matched rows. mu enters both, but on the censored branch only
# through the sum mu + lambdatwo: those rows constrain the total, not the split
# between shared and source-2-only. Separating mu from lambdatwo therefore
# rests on the matched rows as well, and a design with few of them will learn
# the congruence f weakly however many censored rows it has. This is sharper
# than the corresponding statement for binegbin_cens, where the two
# dispersions carry some of the same information, and it is visible directly in
# the closed form above.
#
# EXACT RELATIONSHIP TO THE OTHER FAMILIES. On the matched branch this family
# is bipois exactly -- the R reference below delegates to bipois_lpmf_r rather
# than restating the sum, so the equivalence is structural on the R side and
# tested on the Stan side, where the recurrence is written out separately. As
# shapes, shapexone and shapextwo -> Inf, binegbin_cens -> bipois_cens on
# both branches. Both identities are pinned in test-bipois_cens.R.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Censoring-aware joint bivariate-Poisson family for brms
#'
#' @description
#' Censoring-aware extension of [bipois()]. Models the same trivariate-reduction
#' bivariate-Poisson pair `(y1, y2)` -- `y1 = N_shared + N1`,
#' `y2 = N_shared + N2`, with `N_shared ~ Poisson(mu)`,
#' `N1 ~ Poisson(lambdaone)`, `N2 ~ Poisson(lambdatwo)` mutually independent
#' given their rates -- but allows the first margin (`y1`) to be UNOBSERVED on
#' some rows. Each row carries two supplementary integers via `vint()`: `y2`
#' (the always-observed second count) and `y1_obs` (a 0/1 flag marking whether
#' `y1` was observed for that row).
#'
#' Stands to [bipois()] exactly as [binegbin_cens()] stands to [binegbin()],
#' and is the equidispersed special case of [binegbin_cens()].
#'
#' On `y1_obs == 1` (matched) rows the likelihood is the full joint [bipois()]
#' lpmf on `(y1, y2)`. On `y1_obs == 0` (`y2`-only) rows it is the
#' `y1`-integrated marginal of that same joint, which for Poisson components is
#' available in closed form: a sum of independent Poissons is Poisson, so
#' `y2 ~ Poisson(mu + lambdatwo)` exactly. [binegbin_cens()] must evaluate the
#' corresponding convolution as a sum; this family does not.
#'
#' Three dpars, the same as [bipois()]: `mu` (shared rate), `lambdaone` and
#' `lambdatwo` (the two source-specific rates), all `link = "log"`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2, y1_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaone ~ lamx + methd),
#'        nlf(lambdatwo ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1, nl = TRUE),
#'     family   = bipois_cens(),
#'     stanvars = bipois_cens_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **When to use this rather than [binegbin_cens()].** This family fixes each
#' latent component's variance equal to its mean. Where the counts are genuinely
#' overdispersed relative to that, [binegbin_cens()] is the correct model and
#' this one will understate the marginal variances. Where they are not,
#' [binegbin_cens()] can only represent the fit by driving its dispersions to
#' their Poisson limit (`shape` \eqn{\to\infty}, equivalently `kappa`
#' \eqn{\to 0}), a boundary at which sampling degrades; fitting the
#' equidispersed family directly avoids it. Compare the two with `loo()`.
#'
#' **What each rate is identified from.** `lambdatwo` appears on both branches,
#' so every row informs it. `lambdaone` appears only on the matched branch and is
#' identified solely by the matched rows. `mu` appears on both, but the censored
#' branch sees it only through the sum `mu + lambdatwo` -- those rows constrain
#' the total rate of the observed margin, not how it divides between the shared
#' and source-2-only components. Separating `mu` from `lambdatwo` -- and so
#' estimating the congruence \eqn{f} -- therefore also rests on the matched rows.
#' With few of them, `mu` and `lambdatwo` trade off along their sum and the
#' prior does correspondingly more of the work, however many censored rows the
#' design contains.
#'
#' **Two `vint()` arguments, in declared order.** brms appends `vint()` integers
#' to the generated lpmf call in the order they are listed in the formula's
#' `vint()` term, matching the `vars` declared here
#' (`c("vint1[n]", "vint2[n]")`): so `vint(y2, y1_obs)` binds `vint1 = y2` and
#' `vint2 = y1_obs`. brms generates `target += bipois_cens_lpmf(Y[n] | mu[n],
#' lambdaone[n], lambdatwo[n], vint1[n], vint2[n])` -- dpars in the order
#' declared here, then the two vint args. `bipois_cens_stan_funs` declares
#' `bipois_cens_lpmf` with exactly this signature; reordering the dpars or the
#' two `vint()` terms without matching the Stan signature silently swaps which
#' rate governs which component or which integer is the branch flag.
#'
#' On `y1_obs == 0` rows the response column `Y` is not read by the likelihood,
#' so any integer placeholder there is inert. Supply one rather than `NA`, which
#' brms rejects before the family is reached.
#'
#' **Forced `mu` naming, and the second count via `vint()`.** Identical
#' conventions to [bipois()] -- `mu` is brms's mandatory dpar name, here bound to
#' the shared component's rate, not a mean of either response; `y2` (and
#' `y1_obs`) travel as supplementary integer data through `vint()` because
#' `custom_family()` declares a single response column. See [bipois()] for the
#' full explanation, including why the rates are spelled `lambdaone`/`lambdatwo`
#' in code but written \eqn{\lambda_1}{lambda_1}/\eqn{\lambda_2}{lambda_2} in
#' the documentation.
#'
#' @return A brms custom_family object.
#' @export
bipois_cens <- function() {
  brms::custom_family(
    name  = "bipois_cens",
    dpars = c("mu", "lambdaone", "lambdatwo"),  # mu = lambda_shared -- see bipois()
    links = c("log", "log", "log"),
    lb    = c(0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "vint2[n]")  # vint1 = y2, vint2 = y1_obs
  )
}

#' @rdname bipois_cens
#' @export
bipois_cens_stanvars <- function() {
  brms::stanvar(block = "functions", scode = bipois_cens_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Branching lpmf. The y1_obs == 1 branch is the bipois_lpmf body verbatim --
# the k = 0 term computed directly, each subsequent term obtained from the
# previous by a ratio update in log space, accumulated via log_sum_exp; see
# bipois.R for the derivation of that ratio and its correspondence to the
# cited source. It is written out here rather than calling bipois_lpmf so that
# bipois_cens_stanvars() is self-contained: a model supplying only this
# family's stanvars must still compile.
#
# The y1_obs == 0 branch is the closed-form marginal, poisson_lpmf(y2 | mu +
# lambdatwo). No sum, so nothing to cap or truncate; this is the branch that
# has no counterpart in binegbin_cens, where the same quantity is a
# convolution over k = 0..y2.
bipois_cens_stan_funs <- "
  real bipois_cens_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
                         int y2, int y1_obs) {
    if (y1_obs == 1) {
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
    } else {
      return poisson_lpmf(y2 | mu + lambdatwo);
    }
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Branch-selecting evaluation, vectorised over all arguments (recycled to a
# common length), used to power log_lik_bipois_cens() and to validate the Stan
# lpmf. Internal reference only, mirroring binegbin_cens_lpmf_r's role.
#
# The matched branch DELEGATES to bipois_lpmf_r rather than restating the sum,
# so "the matched branch is bipois" holds by construction on the R side. The
# Stan side writes the recurrence out separately (see above), which is where the
# equivalence is a claim rather than a definition, and that is where the test
# checks it.
#
# The censored branch uses the closed form, as Stan does. The independent route
# -- brute-force convolution over k = 0..y2 -- lives in test-bipois_cens.R
# rather than here, so log_lik() does not pay for a sum that has an exact
# one-line answer. Keeping it in the test is what makes the closed form checked
# rather than assumed.
bipois_cens_lpmf_r <- function(y1, y2, y1_obs, mu, lambdaone, lambdatwo) {
  n <- max(length(y1), length(y2), length(y1_obs), length(mu),
           length(lambdaone), length(lambdatwo))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  y1_obs    <- rep_len(y1_obs, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  out <- numeric(n)
  matched  <- which(y1_obs == 1)
  censored <- which(y1_obs != 1)

  if (length(matched)) {
    out[matched] <- bipois_lpmf_r(
      y1[matched], y2[matched],
      mu[matched], lambdaone[matched], lambdatwo[matched]
    )
  }
  if (length(censored)) {
    out[censored] <- stats::dpois(
      y2[censored], mu[censored] + lambdatwo[censored], log = TRUE
    )
  }
  out
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname bipois_cens
#' @export
#' @keywords internal
log_lik_bipois_cens <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared -- see bipois()
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y1     <- prep$data$Y[i]
  y2     <- prep$data$vint1[i]
  y1_obs <- prep$data$vint2[i]
  bipois_cens_lpmf_r(y1, y2, y1_obs, mu, lambdaone, lambdatwo)
}

#' @rdname bipois_cens
#' @export
#' @keywords internal
posterior_predict_bipois_cens <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)   # lambda_shared -- see bipois()
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y2 <- prep$data$vint1[i]
  # y1_obs is deliberately IGNORED here: every row gets a y1 draw conditional on
  # its observed y2, matched and y2-only alike, which is what imputing the
  # unobserved margin across every row requires. Same convention as
  # posterior_predict_binegbin_cens().
  #
  # Unlike the negative-binomial families, the conditional split is exactly
  # Binomial: conditioning a sum of independent Poissons on its total gives a
  # Binomial split, so N_shared | y2 ~ Binomial(y2, mu / (mu + lambdatwo)) with
  # no discrete conditional to sample. Identical construction to
  # posterior_predict_bipois().
  p_shared <- mu / (mu + lambdatwo)
  n_shared <- stats::rbinom(length(mu), size = y2, prob = p_shared)
  n1       <- stats::rpois(length(mu), lambdaone)
  n_shared + n1
}

#' @rdname bipois_cens
#' @export
#' @keywords internal
posterior_epred_bipois_cens <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")   # lambda_shared -- see bipois()
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + E[N1] = y2 * mu/(mu + lambdatwo) + lambdaone,
  # exact, being the mean of the Binomial split used in
  # posterior_predict_bipois_cens() above. Returned for every row including the
  # censored ones, matching that function's convention: the conditional
  # expectation of the unobserved margin is defined there too, and is the
  # quantity a user imputing y1 wants.
  y2_mat   <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}
