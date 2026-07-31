# ==========================================================================
# binegbin_joint: censoring-aware bivariate Negative-Binomial
#
# The censoring-aware extension of binegbin (see binegbin.R). Same generative
# model -- y1 = N_shared + N1, y2 = N_shared + N2 with
#
#   N_shared ~ NB2(mu,        shapes)   (shared component; drives correlation)
#   N1       ~ NB2(lambdaone, shapex)   (source-1-only excess)
#   N2       ~ NB2(lambdatwo, shapex)   (source-2-only excess)
#
# and the same five dpars and links -- but each row now carries a second
# supplementary integer, an observation flag y1_obs, alongside y2:
#
#   y1_obs == 1 (matched set):  full joint binegbin lpmf on (y1, y2) --
#       identical to binegbin_lpmf, byte-for-byte.
#   y1_obs == 0 (y2-only set):  the y2 MARGINAL of the SAME bivariate model --
#       P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex),
#       i.e. the joint with the y1 (N1) term integrated out over all y1.
#
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS. The unmatched (y2-only)
# rows never observe y1, so a plain binegbin fit could only use the matched
# rows. But the y2-only rows still carry information about the SHARED
# structure (mu, shapes, lambdatwo) and the vessel/trip random effects: their
# y2 is a draw from the same bivariate model, merely with its y1 margin
# unobserved. Integrating y1 out (rather than dropping those rows, or --
# worse -- giving them their own single-dispersion neg_binomial_2 on y2,
# which is a DIFFERENT model inconsistent with the matched decomposition)
# lets one brm() call pool all rows under one coherent likelihood. lambdaone
# (the source-1-only excess rate) and the directional bias between the two
# excess rates are then identified ONLY by the matched rows; the y2-only rows
# sharpen mu / shapes / lambdatwo and the shared vessel+trip structure. This is
# the standard partially-observed-margin construction (a censored/"missing at
# random on the y1 side" likelihood), not a heuristic.
#
# EXACT RELATIONSHIP TO binegbin. On the matched branch this family's lpmf is
# the binegbin lpmf, exactly (same sum, same neg_binomial_2 terms). The
# y2-only branch is the y1-integrated marginal of that same joint. Two
# consequences the package tests pin down:
#   * sum over y1 of the matched-branch lpmf == the y2-only-branch lpmf
#     (marginal identity), and
#   * binegbin_joint_lpmf(y1_obs == 1) == binegbin_lpmf on identical inputs
#     (equivalence).
# The second identity is what licenses reading a binegbin fit's cached draws
# as the y1_obs == 1 slice of a binegbin_joint fit and vice versa -- the two
# families cannot silently drift apart without a test going red.
#
# Validated by the marginal identity, the binegbin equivalence, an
# expose_functions grid cross-check of the Stan lpmf against an independent R
# brute-force reference (~1e-14), and the conditional-prediction identity
# (posterior_predict draws == joint / marginal). See
# tests/testthat/test-binegbin_joint.R.

# --------------------------------------------------------------------------
# brms custom family
# --------------------------------------------------------------------------

#' Censoring-aware joint bivariate-Negative-Binomial family for brms
#'
#' @description
#' Censoring-aware extension of [binegbin()]. Models the same trivariate-
#' reduction bivariate Negative-Binomial pair `(y1, y2)` -- `y1 =
#' N_shared + N1`, `y2 = N_shared + N2`, with `N_shared ~ NB2(mu,
#' shapes)`, `N1 ~ NB2(lambdaone, shapex)`, `N2 ~ NB2(lambdatwo, shapex)`
#' mutually independent given their rates -- but allows the first margin (`y1`)
#' to be UNOBSERVED on some rows. Each row carries two supplementary integers
#' via `vint()`: `y2` (the always-observed second count) and `y1_obs` (a
#' 0/1 flag marking whether `y1` was observed for that row).
#'
#' On `y1_obs == 1` (matched) rows the likelihood is the full joint
#' [binegbin()] lpmf on `(y1, y2)`. On `y1_obs == 0` (`y2`-only) rows it is
#' the `y1`-integrated marginal of that same joint,
#' `P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex)` --
#' NOT a separate single-dispersion `neg_binomial_2` on `y2`, which would be
#' a different model inconsistent with the matched decomposition. This lets
#' one `brm()` call pool matched and `y2`-only rows under one coherent
#' likelihood: `lambdaone` and the between-source bias are identified only by
#' the matched rows, while the `y2`-only rows sharpen `mu`, `shapes`,
#' `lambdatwo`, and the shared vessel/trip random-effect structure.
#'
#' Five dpars, identical to [binegbin()]: the three rates (`mu` = shared rate,
#' `lambdaone`/`lambdatwo` = the two source-specific rates) plus two
#' dispersions -- `shapes` for the shared component and `shapex` shared across
#' the two excess components. All five use `link = "log"` (see [binegbin()]).
#' To share the excess level across the two rates and split them by a
#' directional bias, supply them through non-linear formulas *without* an
#' explicit `exp()` -- the log link applies it, so
#' `nlf(lambdaone ~ lamx + methd)` gives `lambdaone = exp(lamx + methd)`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2, y1_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaone ~ lamx + methd),
#'        nlf(lambdatwo ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1,
#'        shapes ~ 1, shapex ~ 1, nl = TRUE),
#'     family   = binegbin_joint(),
#'     stanvars = binegbin_joint_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Two `vint()` arguments, in declared order.** brms appends `vint()`
#' integers to the generated lpmf call in the order they are listed in the
#' formula's `vint()` term, matching the `vars` declared here
#' (`c("vint1[n]", "vint2[n]")`): so `vint(y2, y1_obs)` binds `vint1 = y2`
#' and `vint2 = y1_obs`. brms generates `target += binegbin_joint_lpmf(Y[n] |
#' mu[n], lambdaone[n], lambdatwo[n], shapes[n], shapex[n], vint1[n], vint2[n])`
#' -- dpars in the order declared here, then the two vint args.
#' `binegbin_joint_stan_funs` declares `binegbin_joint_lpmf` with exactly this
#' signature; reordering the dpars or the two `vint()` terms without matching
#' the Stan signature silently swaps which quantity governs which component or
#' which integer is the branch flag.
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
#' **Relationship to [binegbin()].** On `y1_obs == 1` rows this family's lpmf
#' equals the [binegbin()] lpmf exactly (same marginalisation sum). The
#' `y1_obs == 0` branch is the `y1`-integrated marginal of that same
#' bivariate model. The package tests pin both identities (marginal identity;
#' binegbin equivalence).
#'
#' @return A brms custom_family object.
#' @export
binegbin_joint <- function() {
  brms::custom_family(
    name  = "binegbin_joint",
    dpars = c("mu", "lambdaone", "lambdatwo", "shapes", "shapex"),
    links = c("log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0),
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

# Branching marginalisation sum. The y1_obs == 1 branch is byte-for-byte the
# binegbin_lpmf body (direct log_sum_exp over k = 0..min(y1, y2); see
# binegbin.R for why a direct sum rather than a recurrence). The y1_obs == 0
# branch drops the y1 (N1) term and sums over k = 0..y2, i.e. the same
# joint with the first margin integrated out. neg_binomial_2_lpmf(0 | m, phi)
# is well-defined, so zero counts and the k = 0 term need no special casing.
binegbin_joint_stan_funs <- "
  real binegbin_joint_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
                           real shapes, real shapex, int y2, int y1_obs) {
    if (y1_obs == 1) {
      int m = min(y1, y2);
      vector[m + 1] lp;
      for (k in 0:m) {
        lp[k + 1] = neg_binomial_2_lpmf(k      | mu,        shapes)
                  + neg_binomial_2_lpmf(y1 - k | lambdaone, shapex)
                  + neg_binomial_2_lpmf(y2 - k | lambdatwo, shapex);
      }
      return log_sum_exp(lp);
    } else {
      vector[y2 + 1] lp;
      for (k in 0:y2) {
        lp[k + 1] = neg_binomial_2_lpmf(k      | mu,        shapes)
                  + neg_binomial_2_lpmf(y2 - k | lambdatwo, shapex);
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
binegbin_joint_lpmf_r <- function(y1, y2, y1_obs, mu, lambdaone, lambdatwo,
                                  shapes, shapex) {
  n <- max(length(y1), length(y2), length(y1_obs), length(mu),
           length(lambdaone), length(lambdatwo), length(shapes), length(shapex))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  y1_obs    <- rep_len(y1_obs, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)
  shapes    <- rep_len(shapes, n)
  shapex    <- rep_len(shapex, n)

  vapply(seq_len(n), function(i) {
    if (y1_obs[i] == 1) {
      m <- min(y1[i], y2[i])
      k <- 0:m
      log_terms <- stats::dnbinom(k,      size = shapes[i], mu = mu[i],        log = TRUE) +
        stats::dnbinom(y1[i] - k,         size = shapex[i], mu = lambdaone[i], log = TRUE) +
        stats::dnbinom(y2[i] - k,         size = shapex[i], mu = lambdatwo[i], log = TRUE)
    } else {
      k <- 0:y2[i]
      log_terms <- stats::dnbinom(k,      size = shapes[i], mu = mu[i],        log = TRUE) +
        stats::dnbinom(y2[i] - k,         size = shapex[i], mu = lambdatwo[i], log = TRUE)
    }
    mx <- max(log_terms)
    mx + log(sum(exp(log_terms - mx)))
  }, numeric(1))
}

# --------------------------------------------------------------------------
# brms interface functions -- found by name convention, must be exported
# --------------------------------------------------------------------------

#' @rdname binegbin_joint
#' @export
#' @keywords internal
log_lik_binegbin_joint <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapex    <- brms::get_dpar(prep, "shapex", i = i)
  y1     <- prep$data$Y[i]
  y2     <- prep$data$vint1[i]
  y1_obs <- prep$data$vint2[i]
  binegbin_joint_lpmf_r(y1, y2, y1_obs, mu, lambdaone, lambdatwo, shapes, shapex)
}

#' @rdname binegbin_joint
#' @export
#' @keywords internal
posterior_predict_binegbin_joint <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapex    <- brms::get_dpar(prep, "shapex", i = i)
  y2 <- prep$data$vint1[i]
  # y1_obs is deliberately IGNORED here: every set gets a y1 draw conditional
  # on its observed y2, matched and y2-only alike, which is what imputing the
  # unobserved margin across every row requires. The conditional
  # split N_shared | y2 is NOT Binomial (a NegBin sum condition is not
  # Binomial); it is P(N_shared = k | y2) proportional to
  # NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapex) over k = 0..y2.
  # Sample that discrete conditional, then add a fresh N1 ~ NB2(lambdaone,
  # shapex). Identical construction to posterior_predict_binegbin.
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
