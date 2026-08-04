# ==========================================================================
# binegbin_joint: censoring-aware bivariate Negative-Binomial
#
# The censoring-aware extension of binegbin (see binegbin.R). Same generative
# model -- y1 = N_shared + N1, y2 = N_shared + N2 with
#
#   N_shared ~ NB2(mu,        shapes)      (shared component; drives correlation)
#   N1       ~ NB2(lambdaone, shapexone)   (source-1-only excess)
#   N2       ~ NB2(lambdatwo, shapextwo)   (source-2-only excess)
#
# and six dpars, all log-linked -- but each row now carries a second
# supplementary integer, an observation flag y1_obs, alongside y2:
#
#   y1_obs == 1 (matched set):  full joint binegbin lpmf on (y1, y2).
#   y1_obs == 0 (y2-only set):  the y2 MARGINAL of the SAME bivariate model --
#       P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo),
#       i.e. the joint with the y1 (N1) term integrated out over all y1.
#
# WHY THE TWO EXCESS DISPERSIONS ARE SEPARATE. Up to 0.7.0 a single dpar
# `shapex` governed both excess components, imposing shapexone == shapextwo.
# That constraint is a modelling choice, not a property of the construction:
# the two sources are different instruments and there is no reason their
# source-only excess must be equally overdispersed. 0.8.0 therefore frees the
# two, and the symmetric model becomes a FORMULA constraint rather than a
# separate family -- supply both through one non-linear parameter,
#   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1
# and the fit is term-for-term the pre-0.8.0 five-dpar model. See
# migration/family-unification.md.
#
# Note the asymmetry in what the two dispersions can be identified from.
# shapextwo appears on BOTH branches (it governs the always-observed margin),
# so the y2-only rows sharpen it. shapexone appears only on the matched
# branch, so it is identified SOLELY by the matched rows -- the same way
# lambdaone is. A design with few matched rows will therefore learn shapexone
# far less sharply than shapextwo, and the prior on it does correspondingly
# more of the work.
#
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS. The unmatched (y2-only)
# rows never observe y1, so a plain binegbin fit could only use the matched
# rows. But the y2-only rows still carry information about the SHARED
# structure (mu, shapes, lambdatwo, shapextwo) and the vessel/trip random
# effects: their y2 is a draw from the same bivariate model, merely with its
# y1 margin unobserved. Integrating y1 out (rather than dropping those rows,
# or -- worse -- giving them their own single-dispersion neg_binomial_2 on y2,
# which is a DIFFERENT model inconsistent with the matched decomposition)
# lets one brm() call pool all rows under one coherent likelihood. This is the
# standard partially-observed-margin construction (a censored/"missing at
# random on the y1 side" likelihood), not a heuristic.
#
# EXACT RELATIONSHIP TO binegbin. binegbin carries a single excess dispersion.
# On the matched branch this family's lpmf reduces to the binegbin lpmf
# exactly when shapexone == shapextwo. The y2-only branch is the
# y1-integrated marginal of that same joint. Three consequences the package
# tests pin down:
#   * sum over y1 of the matched-branch lpmf == the y2-only-branch lpmf
#     (marginal identity),
#   * binegbin_joint_lpmf(y1_obs == 1, shapexone == shapextwo) ==
#     binegbin_lpmf on identical inputs (equivalence), and
#   * the six-dpar lpmf evaluated with shapexone == shapextwo reproduces the
#     pre-0.8.0 five-dpar lpmf to machine precision (non-regression).
#
# Validated by the marginal identity, the binegbin equivalence, an
# expose_functions grid cross-check of the Stan lpmf against an independent R
# brute-force reference (~1e-14), the conditional-prediction identity
# (posterior_predict draws == joint / marginal), and a simulation-recovery
# check with shapexone far from shapextwo. See
# tests/testthat/test-binegbin_joint.R and test-binegbin_joint_asym.R.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Censoring-aware joint bivariate-Negative-Binomial family for brms
#'
#' @description
#' Censoring-aware extension of [binegbin()]. Models the same trivariate-
#' reduction bivariate Negative-Binomial pair `(y1, y2)` -- `y1 =
#' N_shared + N1`, `y2 = N_shared + N2`, with `N_shared ~ NB2(mu,
#' shapes)`, `N1 ~ NB2(lambdaone, shapexone)`, `N2 ~ NB2(lambdatwo,
#' shapextwo)` mutually independent given their rates -- but allows the first
#' margin (`y1`) to be UNOBSERVED on some rows. Each row carries two
#' supplementary integers via `vint()`: `y2` (the always-observed second count)
#' and `y1_obs` (a 0/1 flag marking whether `y1` was observed for that row).
#'
#' On `y1_obs == 1` (matched) rows the likelihood is the full joint
#' [binegbin()] lpmf on `(y1, y2)`. On `y1_obs == 0` (`y2`-only) rows it is
#' the `y1`-integrated marginal of that same joint,
#' `P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo)` --
#' NOT a separate single-dispersion `neg_binomial_2` on `y2`, which would be
#' a different model inconsistent with the matched decomposition. This lets
#' one `brm()` call pool matched and `y2`-only rows under one coherent
#' likelihood: `lambdaone`, the between-source bias and `shapexone` are
#' identified only by the matched rows, while the `y2`-only rows sharpen `mu`,
#' `shapes`, `lambdatwo`, `shapextwo`, and the shared vessel/trip
#' random-effect structure.
#'
#' Six dpars: the three rates (`mu` = shared rate, `lambdaone`/`lambdatwo` =
#' the two source-specific rates) plus three dispersions -- `shapes` for the
#' shared component and `shapexone`/`shapextwo` for the two source-specific
#' excess components. All six use `link = "log"` (see [binegbin()]). To share a
#' level across a pair of dpars and split them by a directional bias, supply
#' them through non-linear formulas *without* an explicit `exp()` -- the log
#' link applies it, so `nlf(lambdaone ~ lamx + methd)` gives
#' `lambdaone = exp(lamx + methd)`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2, y1_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaone ~ lamx + methd),
#'        nlf(lambdatwo ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1,
#'        shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE),
#'     family   = binegbin_joint(),
#'     stanvars = binegbin_joint_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **The symmetric model is a formula constraint.** Before 0.8.0 this family
#' carried a single excess dispersion `shapex` shared by both margins. That
#' model is the constraint `shapexone == shapextwo`, obtained by routing both
#' through one non-linear parameter:
#'
#' ```
#' bf(y1 | vint(y2, y1_obs) ~ 1,
#'    mu ~ 1 + (1 | vessel),
#'    nlf(shapexone ~ shapexx),
#'    nlf(shapextwo ~ shapexx),
#'    shapexx ~ 1, ..., nl = TRUE)
#' ```
#'
#' The resulting likelihood is term-for-term the pre-0.8.0 five-dpar one (a
#' package test pins this). Stored five-dpar fits keep working unchanged: their
#' single `shapex` resolves to both `shapexone` and `shapextwo` when
#' post-processed. See `migration/family-unification.md` for the migration detail.
#'
#' **What each dispersion is identified from.** `shapextwo` governs the
#' always-observed margin and so appears on both branches; the `y2`-only rows
#' inform it. `shapexone` appears only on the matched branch and is identified
#' solely by the matched rows, as `lambdaone` is. With few matched rows the
#' posterior for `shapexone` is correspondingly dominated by its prior.
#'
#' **Two `vint()` arguments, in declared order.** brms appends `vint()`
#' integers to the generated lpmf call in the order they are listed in the
#' formula's `vint()` term, matching the `vars` declared here
#' (`c("vint1[n]", "vint2[n]")`): so `vint(y2, y1_obs)` binds `vint1 = y2`
#' and `vint2 = y1_obs`. brms generates `target += binegbin_joint_lpmf(Y[n] |
#' mu[n], lambdaone[n], lambdatwo[n], shapes[n], shapexone[n], shapextwo[n],
#' vint1[n], vint2[n])` -- dpars in the order declared here, then the two vint
#' args. `binegbin_joint_stan_funs` declares `binegbin_joint_lpmf` with exactly
#' this signature; reordering the dpars or the two `vint()` terms without
#' matching the Stan signature silently swaps which quantity governs which
#' component or which integer is the branch flag.
#'
#' **Forced `mu` naming, and the second count via `vint()`.** Identical
#' conventions to [binegbin()]/[bipois()] -- `mu` is brms's mandatory dpar
#' name, here bound to the shared component's rate, not a mean of either
#' response; `y2` (and `y1_obs`) travel as supplementary integer data
#' through `vint()` because `custom_family()` declares a single response
#' column. See [bipois()] for the full explanation, including why the rates are
#' spelled `lambdaone`/`lambdatwo` in code but written
#' \eqn{\lambda_1}{lambda_1}/\eqn{\lambda_2}{lambda_2} in the documentation.
#'
#' **Relationship to [binegbin()].** On `y1_obs == 1` rows with
#' `shapexone == shapextwo` this family's lpmf equals the [binegbin()] lpmf
#' exactly (same marginalisation sum). The `y1_obs == 0` branch is the
#' `y1`-integrated marginal of that same bivariate model. The package tests pin
#' both identities (marginal identity; binegbin equivalence).
#'
#' @return A brms custom_family object.
#' @export
binegbin_joint <- function() {
  brms::custom_family(
    name  = "binegbin_joint",
    dpars = c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo"),
    links = c("log", "log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "vint2[n]")  # vint1 = y2, vint2 = y1_obs
  )
}

#' @rdname binegbin_joint
#' @export
binegbin_joint_stanvars <- function() {
  brms::stanvar(block = "functions", scode = binegbin_joint_stan_funs)
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Branching marginalisation sum. The y1_obs == 1 branch is the binegbin_lpmf
# body with the two excess terms carrying their own dispersions (direct
# log_sum_exp over k = 0..min(y1, y2); see binegbin.R for why a direct sum
# rather than a recurrence). The y1_obs == 0 branch drops the y1 (N1) term --
# and with it shapexone -- and sums over k = 0..y2, i.e. the same joint with
# the first margin integrated out. neg_binomial_2_lpmf(0 | m, phi) is
# well-defined, so zero counts and the k = 0 term need no special casing.
binegbin_joint_stan_funs <- "
  real binegbin_joint_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
                           real shapes, real shapexone, real shapextwo,
                           int y2, int y1_obs) {
    if (y1_obs == 1) {
      int m = min(y1, y2);
      vector[m + 1] lp;
      for (k in 0:m) {
        lp[k + 1] = neg_binomial_2_lpmf(k      | mu,        shapes)
                  + neg_binomial_2_lpmf(y1 - k | lambdaone, shapexone)
                  + neg_binomial_2_lpmf(y2 - k | lambdatwo, shapextwo);
      }
      return log_sum_exp(lp);
    } else {
      vector[y2 + 1] lp;
      for (k in 0:y2) {
        lp[k + 1] = neg_binomial_2_lpmf(k      | mu,        shapes)
                  + neg_binomial_2_lpmf(y2 - k | lambdatwo, shapextwo);
      }
      return log_sum_exp(lp);
    }
  }
"

# --------------------------------------------------------------------------
# R-side reference implementation
# --------------------------------------------------------------------------

# Independent brute-force evaluation of the same branching sum via R's
# dnbinom, used to validate the Stan lpmf and to power
# log_lik_binegbin_joint()/posterior_predict_binegbin_joint() post-hoc.
# Internal reference only, mirroring binegbin_lpmf_r's role. Vectorised over
# all arguments (recycled to common length); y1_obs selects the branch
# per-row.
#
# `shapextwo` defaults to `shapexone`, so an eight-argument call is the
# symmetric (pre-0.8.0) model and needs no rewriting.
binegbin_joint_lpmf_r <- function(y1, y2, y1_obs, mu, lambdaone, lambdatwo,
                                  shapes, shapexone, shapextwo = shapexone) {
  n <- max(length(y1), length(y2), length(y1_obs), length(mu),
           length(lambdaone), length(lambdatwo), length(shapes),
           length(shapexone), length(shapextwo))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  y1_obs    <- rep_len(y1_obs, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)
  shapes    <- rep_len(shapes, n)
  shapexone <- rep_len(shapexone, n)
  shapextwo <- rep_len(shapextwo, n)

  vapply(seq_len(n), function(i) {
    if (y1_obs[i] == 1) {
      m <- min(y1[i], y2[i])
      k <- 0:m
      log_terms <- stats::dnbinom(k, size = shapes[i],    mu = mu[i],        log = TRUE) +
        stats::dnbinom(y1[i] - k,   size = shapexone[i],  mu = lambdaone[i], log = TRUE) +
        stats::dnbinom(y2[i] - k,   size = shapextwo[i],  mu = lambdatwo[i], log = TRUE)
    } else {
      k <- 0:y2[i]
      log_terms <- stats::dnbinom(k, size = shapes[i],    mu = mu[i],        log = TRUE) +
        stats::dnbinom(y2[i] - k,   size = shapextwo[i],  mu = lambdatwo[i], log = TRUE)
    }
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

# Accepted spellings for each excess dispersion, most recent first. A stored
# fit is post-processed under the names ITS OWN family object declares, so all
# three eras resolve: 0.8.0+ (shapexone/shapextwo), the project-local
# asymmetric family (shapexem/shapexlb), and <= 0.7.0 symmetric fits, whose
# single `shapex` correctly serves both margins. See .get_dpar_any() in
# utils.R.
.SHAPEXONE_NAMES <- c("shapexone", "shapexem", "shapex")
.SHAPEXTWO_NAMES <- c("shapextwo", "shapexlb", "shapex")

#' @rdname binegbin_joint
#' @export
#' @keywords internal
log_lik_binegbin_joint <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y1     <- prep$data$Y[i]
  y2     <- prep$data$vint1[i]
  y1_obs <- prep$data$vint2[i]
  binegbin_joint_lpmf_r(y1, y2, y1_obs, mu, lambdaone, lambdatwo,
                        shapes, shapexone, shapextwo)
}

#' @rdname binegbin_joint
#' @export
#' @keywords internal
posterior_predict_binegbin_joint <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y2 <- prep$data$vint1[i]
  # y1_obs is deliberately IGNORED here: every set gets a y1 draw conditional
  # on its observed y2, matched and y2-only alike, which is what imputing the
  # unobserved margin across every row requires. The conditional
  # split N_shared | y2 is NOT Binomial (a NegBin sum condition is not
  # Binomial); it is P(N_shared = k | y2) proportional to
  # NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo) over k = 0..y2 --
  # the y2 margin, hence shapextwo. Sample that discrete conditional, then add
  # a fresh N1 ~ NB2(lambdaone, shapexone). Identical construction to
  # posterior_predict_binegbin.
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

#' @rdname binegbin_joint
#' @export
#' @keywords internal
posterior_epred_binegbin_joint <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  shapes    <- brms::get_dpar(prep, "shapes")
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES)
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + lambdaone, the expectation of exactly what
  # posterior_predict_binegbin_joint() simulates. `shapextwo` and not
  # `shapexone`: the quantity conditioned on is y2, so it is the SECOND
  # margin's excess dispersion that enters the conditional weights. See
  # .e_shared_given_y2_nb() in utils.R.
  #
  # y1_obs is ignored, as it is in posterior_predict_binegbin_joint(): the
  # conditional expectation of the first margin is defined on censored rows too,
  # and imputing it there is the point of having the family. Returning it for
  # every row keeps epred and posterior_predict comparable row by row.
  #
  # This family had no posterior_epred before 0.9.0, on the view that E[y1] is
  # ambiguous when y1 is unobserved. It is not: E[y1 | y2] is well defined
  # whether or not y1 was recorded, and it is the quantity a user imputing the
  # missing margin wants.
  .e_shared_given_y2_nb(mu, lambdatwo, shapes, shapextwo, y2) + lambdaone
}
