# ==========================================================================
# binegbin: joint bivariate Negative-Binomial via trivariate reduction
#
# The overdispersed sibling of bipois (see bipois.R). Same trivariate-
# reduction construction -- y1 = N_shared + N1, y2 = N_shared + N2,
# the three latent counts mutually independent given their rates -- but each
# latent count is Negative-Binomial rather than Poisson:
#
#   N_shared ~ NB2(mu,        shapes)      (shared component; drives correlation)
#   N1       ~ NB2(lambdaone, shapexone)   (source-1-only excess)
#   N2       ~ NB2(lambdatwo, shapextwo)   (source-2-only excess)
#
# NB2(m, phi) is Stan's neg_binomial_2 / R's dnbinom(size = phi, mu = m):
# mean m, variance m + m^2/phi. shapes is the shared-component dispersion;
# shapexone and shapextwo are the two source-specific excess dispersions.
#
# ONE FAMILY, TWO CONSTRUCTORS. binegbin() is for a fully paired design and
# binegbin_partialobs() for one in which the first count is missing on some
# rows. Both return the SAME custom_family name, so there is one lpmf and one
# set of post-processing methods; they differ only in how the observation flag
# reaches the likelihood:
#
#   binegbin()             vars = c("vint1[n]", "1")
#   binegbin_partialobs()  vars = c("vint1[n]", "vint2[n]")
#
# brms pastes `vars` entries into the generated Stan call verbatim, so the
# plain constructor supplies y1_obs as the literal 1 and its user never has to
# build a flag column, or know one exists. See binegbin_partialobs() for why
# this is preferred over declaring two overloaded Stan functions.
#
# WHAT THE FLAG SELECTS.
#
#   y1_obs == 1 (matched row):  full joint lpmf on (y1, y2).
#   y1_obs == 0 (y2-only row):  the y2 MARGINAL of the SAME bivariate model --
#       P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo),
#       i.e. the joint with the y1 (N1) term integrated out over all y1.
#
# This is not censoring in brms's sense. brms's own cens() addition term means
# a value known to lie in a set; here y1 is not observed at all and the
# likelihood marginalises over its whole support. The families were named
# `_cens` up to 0.9.1 for that reason and renamed at 0.10.0.
#
# WHY BOTH DISPERSIONS ARE FREE IN BOTH CONSTRUCTORS. Up to 0.7.0 a single
# dpar `shapex` governed both excess components, imposing
# shapexone == shapextwo. That constraint is a modelling choice, not a
# property of the construction: the two sources are different instruments and
# there is no reason their source-only excess must be equally overdispersed.
# 0.8.0 freed the two for the partially observed family only; 0.10.0 does so
# for both, which is the right way round. Every row of a fully paired design
# informs both dispersions, so that is the case in which they are EASIEST to
# identify; it was the case that lacked the capability.
#
# The symmetric model is a FORMULA constraint rather than a separate family --
# supply both through one non-linear parameter,
#   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1
# and the fit is term-for-term the pre-0.8.0 five-dpar model. See
# migration/family-unification.md.
#
# Note the asymmetry in what the two dispersions can be identified from under
# PARTIAL observation. shapextwo appears on both branches (it governs the
# always-observed margin), so the y2-only rows sharpen it. shapexone appears
# only on the matched branch, so it is identified SOLELY by the matched rows
# -- the same way lambdaone is. A design with few matched rows will therefore
# learn shapexone far less sharply than shapextwo, and the prior on it does
# correspondingly more of the work. Under full pairing the asymmetry vanishes.
#
# WHY NEGBIN AND NOT AN OLRE ON bipois. The plain-Poisson bipois cannot be
# overdispersed (Var == mean for each latent count), so it underfits the
# real marginal variances badly -- in one motivating dataset by ~10x
# (Var(y1) fitted 16.6 against 179 observed) and Var(d) by ~3.5x. The
# obvious fix -- add a per-set observation-level random effect (OLRE) on the
# excess components -- FAILS synthetic recovery: with one bivariate
# observation per set but three per-set latent deviates (mu-OLRE + two excess
# OLREs), the excess deviates act as residual-absorbers, their population SD
# collapses toward the prior mode, and drawing fresh deviates does NOT
# regenerate the observed spread (recovered excess SD 0.37 vs true 0.85;
# fresh-deviate Var(d) 2.9 vs true 19.2). A conditional posterior-predictive
# check hides this completely -- only a marginal (fresh-deviate) check exposes
# it. binegbin carries the dispersion in SCALAR shapes/shapexone/shapextwo
# instead, estimated from aggregate mean-variance mismatch across sets --
# identifiable, no per-set overfitting, clean marginal PPC, and consistent
# with the review-track NegBin models.
#
# LIKELIHOOD. N_shared is unobserved and marginalised out analytically,
# exactly as in bipois -- the sum structure is identical, only the component
# pmfs change from Poisson to NegBin:
#
#   P(y1=x, y2=y) = sum_{k=0}^{min(x,y)}
#     NB2(k | mu, shapes) NB2(x-k | lambdaone, shapexone)
#                         NB2(y-k | lambdatwo, shapextwo)
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
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS, under partial observation.
# The unmatched (y2-only) rows never observe y1, so a matched-only fit could
# use the matched rows alone. But the y2-only rows still carry information
# about the SHARED structure (mu, shapes, lambdatwo, shapextwo) and the
# vessel/trip random effects: their y2 is a draw from the same bivariate
# model, merely with its y1 margin unobserved. Integrating y1 out (rather than
# dropping those rows, or -- worse -- giving them their own single-dispersion
# neg_binomial_2 on y2, which is a DIFFERENT model inconsistent with the
# matched decomposition) lets one brm() call pool all rows under one coherent
# likelihood. This is the standard partially-observed-margin construction, not
# a heuristic.
#
# Validation -- grid cross-check of the Stan lpmf against the independent R
# brute-force reference to ~1e-14, normalisation to 1 on both branches, the
# moment identities, the marginal identity (summing the matched branch over y1
# gives the y2-only branch), the Poisson-limit reduction to bipois, the
# conditional-prediction identity (posterior_predict draws == joint /
# marginal), and end-to-end parameter recovery with a coverage assessment --
# is in tests/testthat/test-binegbin.R, test-binegbin-partialobs.R and
# test-partialobs-predict.R.

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
#' `N_shared ~ NB2(mu, shapes)`, `N1 ~ NB2(lambdaone, shapexone)`,
#' `N2 ~ NB2(lambdatwo, shapextwo)` mutually independent given their rates.
#' `NB2(m, phi)` has mean `m` and variance `m + m^2/phi` (Stan
#' `neg_binomial_2`; R `dnbinom(size = phi, mu = m)`).
#'
#' Use this when both counts were recorded on every row. If the first count is
#' missing on some rows, use [binegbin_partialobs()], which is the same family
#' with an observation flag -- same `name`, same six dpars, same likelihood,
#' same post-processing.
#'
#' Six dpars: the three rates (`mu` = shared rate, `lambdaone`/`lambdatwo` =
#' the two source-specific rates) plus three dispersions -- `shapes` for the
#' shared component and `shapexone`/`shapextwo` for the two source-specific
#' excess components. All six use `link = "log"`. Supply the excess rates
#' through a non-linear formula without an explicit `exp()` (the log link
#' applies it): `nlf(lambdaone ~ lamx)` gives `lambdaone = exp(lamx)`.
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
#'        shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE),
#'     family   = binegbin(),
#'     stanvars = binegbin_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **Two excess dispersions, and the symmetric special case.** Up to 0.9.1
#' this family carried a single `shapex` shared by both excess components,
#' imposing `shapexone == shapextwo`. That is a modelling choice rather than a
#' property of the construction, and 0.10.0 frees it. To recover the
#' constraint, route both through one non-linear parameter:
#'
#' ```
#' bf(y1 | vint(y2) ~ 1,
#'    mu ~ 1 + (1 | vessel),
#'    nlf(shapexone ~ shapexx),
#'    nlf(shapextwo ~ shapexx),
#'    shapexx ~ 1, ..., nl = TRUE)
#' ```
#'
#' The resulting likelihood is term-for-term the pre-0.10.0 five-dpar one, and
#' a package test pins that. Stored five-dpar fits keep working unchanged:
#' their single `shapex` resolves to both `shapexone` and `shapextwo` when
#' post-processed. Set a prior with `nlpar = "shapexx"` rather than
#' `dpar = "shapex"`.
#'
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
#' shapes[n], shapexone[n], shapextwo[n], vint1[n], 1)` -- dpars in the order
#' declared here, then the two `vars` entries. `binegbin_stan_funs` declares
#' `binegbin_lpmf` with exactly this signature; reordering one without the
#' other silently swaps which rate or dispersion governs which component.
#' The trailing `1` is the observation flag, supplied as a literal because
#' every row of a fully paired design has both counts.
#'
#' @return A brms custom_family object.
#' @seealso [binegbin_partialobs()] for partially observed pairs; [bipois()]
#'   for the equidispersed Poisson counterpart.
#' @export
binegbin <- function() {
  brms::custom_family(
    name  = "binegbin",
    dpars = c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo"),
    links = c("log", "log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "1")  # vint1 = y2; the flag is a literal -- see Details
  )
}

#' @rdname binegbin
#' @export
binegbin_stanvars <- function() {
  brms::stanvar(block = "functions", scode = binegbin_stan_funs)
}

#' Joint bivariate-Negative-Binomial family for partially observed pairs
#'
#' @description
#' [binegbin()] for a design in which the first count is missing on some rows.
#' Same generative model, same `name`, same six dpars, same likelihood and the
#' same post-processing methods -- the only difference is that each row carries
#' a second supplementary integer, an observation flag, through `vint()`:
#'
#' ```
#' bf(y1 | vint(y2, y1_obs) ~ ...)
#' ```
#'
#' `y1_obs` is a 0/1 integer column: `1` where both counts were recorded, `0`
#' where the first was not. `y1` may hold any non-negative integer on those
#' rows -- `0` is the conventional placeholder -- because the likelihood does
#' not read it. Do not use `NA`, which brms drops before fitting, taking the
#' row's observed `y2` with it.
#'
#' **What happens to each kind of row.** A matched row (`y1_obs == 1`) uses the
#' full joint lpmf on `(y1, y2)`. A row whose first count was never recorded
#' (`y1_obs == 0`) contributes the second count's marginal *from the same
#' model*,
#' `P(y2) = sum_k NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo)` --
#' the joint with the `y1` term integrated out over its whole support. It is
#' not dropped, and it is not given a different model: it still informs the
#' shared component (`mu`, `shapes`), the second source's rate and dispersion
#' (`lambdatwo`, `shapextwo`), and any group-level effects.
#'
#' **What you get afterwards.** The fitted model can impute the unobserved
#' first count conditional on the observed second one, which is usually why
#' someone wanted this. `posterior_predict()` and `posterior_epred()` return a
#' `y1` draw and `E[y1 | y2]` for *every* row, matched and unmatched alike --
#' `y1_obs` selects a likelihood branch, not a prediction.
#'
#' **The design consequence, worth knowing before collecting data.** The first
#' source's rate `lambdaone` and excess dispersion `shapexone` appear only on
#' the matched branch, so they are identified by the matched rows *alone*. A
#' design with 20 matched rows in 500 learns them weakly and leans on their
#' priors. `mu`, `shapes`, `lambdatwo` and `shapextwo` appear on both branches
#' and are informed by every row.
#'
#' **This is not censoring in brms's sense.** brms's `cens()` addition term
#' means a value known to lie in a set -- `left`, `right`, `interval`. Here the
#' first count is not observed at all and the likelihood marginalises over its
#' whole support. This family was called `binegbin_cens()` up to 0.9.1; the
#' name was wrong and was changed at 0.10.0. Do not combine this family with
#' `cens()`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2, y1_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaone ~ lamx + methd),
#'        nlf(lambdatwo ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1,
#'        shapes ~ 1, shapexone ~ 1, shapextwo ~ 1, nl = TRUE),
#'     family   = binegbin_partialobs(),
#'     stanvars = binegbin_partialobs_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **One likelihood, two constructors.** This returns the same
#' `custom_family` `name` as [binegbin()], so both land on one
#' `binegbin_lpmf` and one set of `log_lik_binegbin()` /
#' `posterior_predict_binegbin()` / `posterior_epred_binegbin()` methods. The
#' matched branch of the likelihood is therefore not a second copy that can
#' drift from the fully paired one; it is the same code.
#'
#' What differs is `vars`. [binegbin()] declares `c("vint1[n]", "1")` and this
#' declares `c("vint1[n]", "vint2[n]")`. brms pastes those entries into the
#' generated call, so the fully paired model reaches the same Stan function
#' with the flag fixed at `1`. The alternative -- two Stan functions of the
#' same name and different arity -- would need user-defined function
#' overloading, which arrived in Stan 2.29 (February 2022) and would oblige
#' this package to declare a floor on the Stan version. It does not.
#'
#' **Two `vint()` arguments, in declared order.** brms appends `vint()`
#' integers to the generated lpmf call in the order they are listed in the
#' formula's `vint()` term, matching the `vars` declared here: so
#' `vint(y2, y1_obs)` binds `vint1 = y2` and `vint2 = y1_obs`. Reordering the
#' two `vint()` terms without matching the Stan signature silently swaps the
#' second count with the branch flag.
#'
#' **Telling which shape a stored fit used.** `family$name` is `"binegbin"`
#' either way. What distinguishes them is the presence of the second
#' supplementary integer:
#'
#' ```
#' "vint2" %in% names(brms::standata(fit))   # TRUE for a partially observed fit
#' fit$family$vars     # c("vint1[n]", "vint2[n]") or c("vint1[n]", "1")
#' ```
#'
#' @return A brms custom_family object.
#' @seealso [binegbin()] for the fully paired case; [bipois_partialobs()] for
#'   the equidispersed Poisson counterpart.
#' @export
binegbin_partialobs <- function() {
  brms::custom_family(
    name  = "binegbin",
    dpars = c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo"),
    links = c("log", "log", "log", "log", "log", "log"),
    lb    = c(0, 0, 0, 0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "vint2[n]")  # vint1 = y2, vint2 = y1_obs
  )
}

#' @rdname binegbin_partialobs
#' @export
binegbin_partialobs_stanvars <- function() {
  # Deliberately the same object as binegbin_stanvars(): one family name means
  # one Stan function, and both constructors need exactly it. The alias exists
  # so that a call reads consistently, not because there is a second
  # implementation.
  binegbin_stanvars()
}

# --------------------------------------------------------------------------
# Stan function block
# --------------------------------------------------------------------------

# Branching marginalisation sum. The y1_obs == 1 branch is the direct
# log_sum_exp over k = 0..min(y1, y2) (see the file header for why a direct
# sum rather than a recurrence), each term a sum of three neg_binomial_2 log
# densities. The y1_obs == 0 branch drops the y1 (N1) term -- and with it
# shapexone -- and sums over k = 0..y2, i.e. the same joint with the first
# margin integrated out. neg_binomial_2_lpmf(0 | m, phi) is well-defined, so
# zero counts and the k = 0 term need no special casing.
#
# binegbin() reaches this function with y1_obs supplied as the literal 1, so a
# fully paired model always takes the first branch. There is one lpmf rather
# than two because the matched branch of a partially observed model IS the
# fully paired likelihood, and keeping a second copy is how the two drift.
binegbin_stan_funs <- "
  real binegbin_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
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

# Independent brute-force evaluation of the same branching sum via R's dnbinom
# (an independent route from Stan's neg_binomial_2), used to validate the Stan
# lpmf and to power log_lik_binegbin() post-hoc. Internal reference only,
# mirroring bipois_lpmf_r's role. Vectorised over all arguments (recycled to
# common length); y1_obs selects the branch per row.
#
# Two defaults keep the common calls short and are load-bearing for the test
# suite: `shapextwo` defaults to `shapexone`, so a seven-argument call is the
# symmetric (pre-0.8.0) model, and `y1_obs` defaults to 1, so a call that
# names no flag is the matched branch.
binegbin_lpmf_r <- function(y1, y2, mu, lambdaone, lambdatwo, shapes,
                            shapexone, shapextwo = shapexone, y1_obs = 1L) {
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
# asymmetric family (shapexem/shapexlb), and <= 0.9.1 symmetric binegbin fits,
# whose single `shapex` correctly serves both margins. See .get_dpar_any() in
# utils.R.
.SHAPEXONE_NAMES <- c("shapexone", "shapexem", "shapex")
.SHAPEXTWO_NAMES <- c("shapextwo", "shapexlb", "shapex")

#' @rdname binegbin
#' @export
#' @keywords internal
log_lik_binegbin <- function(i, prep) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y1     <- prep$data$Y[i]
  y2     <- prep$data$vint1[i]
  y1_obs <- .y1_obs_at(prep, i)
  binegbin_lpmf_r(y1, y2, mu, lambdaone, lambdatwo,
                  shapes, shapexone, shapextwo, y1_obs)
}

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_predict_binegbin <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y2 <- prep$data$vint1[i]
  # y1_obs is deliberately NOT read here: every row gets a y1 draw conditional
  # on its observed y2, matched and y2-only alike, which is what imputing the
  # unobserved margin across every row requires. The conditional split
  # N_shared | y2 is NOT Binomial (a NegBin sum condition is not Binomial); it
  # is P(N_shared = k | y2) proportional to
  # NB2(k | mu, shapes) NB2(y2 - k | lambdatwo, shapextwo) over k = 0..y2 --
  # the y2 margin, hence shapextwo. Sample that discrete conditional, then add
  # a fresh N1 ~ NB2(lambdaone, shapexone).
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

#' @rdname binegbin
#' @export
#' @keywords internal
posterior_epred_binegbin <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")        # lambda_shared
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  shapes    <- brms::get_dpar(prep, "shapes")
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES)
  y2 <- prep$data$vint1
  # E[y1 | y2] = E[N_shared | y2] + lambdaone, exact, and the expectation of
  # exactly what posterior_predict_binegbin() simulates. A sum of independent
  # negative binomials conditioned on its total is not Binomial, so
  # E[N_shared | y2] has no closed form as clean as bipois's
  # y2 * mu/(mu+lambdatwo); it is the mean of the discrete conditional over
  # k = 0..y2 -- the same weights posterior_predict_binegbin() samples from,
  # summed rather than sampled. See .e_shared_given_y2_nb() in utils.R.
  #
  # `shapextwo` and not `shapexone`: the quantity conditioned on is y2, so it
  # is the SECOND margin's excess dispersion that enters the conditional
  # weights.
  #
  # y1_obs is not read, as it is not in posterior_predict_binegbin(): the
  # conditional expectation of the first margin is defined on unmatched rows
  # too, and imputing it there is the point of having the family. Returning it
  # for every row keeps epred and posterior_predict comparable row by row.
  #
  # Before 0.9.0 this substituted the MARGINAL shared fraction
  # mu/(mu+lambdatwo), which is bipois's answer and exact only in the Poisson
  # limit, so the epred and posterior_predict disagreed by more than Monte
  # Carlo error. They now agree by construction.
  .e_shared_given_y2_nb(mu, lambdatwo, shapes, shapextwo, y2) + lambdaone
}
