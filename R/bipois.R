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
#
# ONE FAMILY, TWO CONSTRUCTORS. bipois() is for a fully paired design and
# bipois_partialobs() for one in which the first count is missing on some
# rows. Both return the SAME custom_family name, so there is one lpmf and one
# set of post-processing methods; they differ only in how the observation flag
# reaches the likelihood:
#
#   bipois()             vars = c("vint1[n]", "1")
#   bipois_partialobs()  vars = c("vint1[n]", "vint2[n]")
#
# brms pastes `vars` entries into the generated Stan call verbatim, so the
# plain constructor supplies y1_obs as the literal 1 and its user never has to
# build a flag column, or know one exists. See binegbin_partialobs() for why
# this is preferred over declaring two overloaded Stan functions.
#
# WHAT THE FLAG SELECTS.
#
#   y1_obs == 1 (matched row):  full joint bipois lpmf on (y1, y2).
#   y1_obs == 0 (y2-only row):  the y2 MARGINAL of the SAME bivariate model.
#
# This is not censoring in brms's sense: brms's own cens() addition term means
# a value known to lie in a set, whereas here y1 is not observed at all and
# the likelihood marginalises over its whole support. The partially observed
# family was called bipois_cens() up to 0.9.1 and renamed at 0.10.0.
#
# THE UNMATCHED BRANCH IS ANALYTIC, AND THAT IS THE POINT. For binegbin the
# y1-integrated marginal is a convolution that must be evaluated as a sum,
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
# consequences beyond the branch being cheap to evaluate:
#
#   * It supplies an ANALYTIC REFERENCE for binegbin's unmatched branch, which
#     the package otherwise lacks. That branch's marginal-identity test checks
#     the NB sum against itself -- summing the matched branch over y1 and
#     comparing to the unmatched branch -- so both sides share any error in
#     the convolution. Taking the NB family to its Poisson limit and comparing
#     against this closed form is an independent check of the same code path.
#     See test-bipois-partialobs.R.
#   * It avoids a boundary. A user whose counts really are equidispersed must
#     otherwise fit binegbin with the dispersions pressed against their
#     Poisson limit (shape -> Inf, i.e. kappa -> 0), which is exactly where
#     sampling degrades. A dedicated family sidesteps that.
#
# WHY A DEDICATED FAMILY AND NOT TWO SEPARATE FITS, under partial observation.
# Unchanged from binegbin, and the argument does not depend on the component
# distribution: the unmatched rows never observe y1, so a matched-only fit
# could use the matched rows alone, yet those rows' y2 is a draw from the same
# bivariate model and still informs mu, lambdatwo and any group-level
# structure. Integrating the unobserved margin out pools every row under one
# coherent likelihood. See binegbin.R for the full statement.
#
# WHAT EACH RATE IS IDENTIFIED FROM, under partial observation. lambdatwo
# enters both branches, so every row sharpens it. lambdaone enters only the
# matched branch and is identified solely by the matched rows. mu enters both,
# but on the unmatched branch only through the sum mu + lambdatwo: those rows
# constrain the total, not the split between shared and source-2-only.
# Separating mu from lambdatwo therefore rests on the matched rows as well, and
# a design with few of them will learn the congruence f weakly however many
# unmatched rows it has. This is sharper than the corresponding statement for
# binegbin, where the two dispersions carry some of the same information, and
# it is visible directly in the closed form above. Under full pairing none of
# this applies.
#
# EXACT RELATIONSHIP TO binegbin. As shapes, shapexone and shapextwo -> Inf,
# binegbin -> bipois on both branches. Pinned in test-bipois-partialobs.R.

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
#' Use this when both counts were recorded on every row. If the first count is
#' missing on some rows, use [bipois_partialobs()], which is the same family
#' with an observation flag -- same `name`, same three dpars, same likelihood,
#' same post-processing.
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
#' per-row data, consistent with the fact that every row this constructor is
#' for comes from the matched (both-observed) subset.
#'
#' **Order of dpars matters for the generated Stan call.** brms generates
#' `target += bipois_lpmf(Y[n] | mu[n], lambdaone[n], lambdatwo[n],
#' vint1[n], 1)` -- dpars in the order declared here, then the two `vars`
#' entries. `bipois_stan_funs` (stanfunctions via `bipois_stanvars()`)
#' declares `bipois_lpmf` with exactly this argument order; changing the order
#' here without changing the Stan signature (or vice versa) silently swaps
#' which rate governs which count. The trailing `1` is the observation flag,
#' supplied as a literal because every row of a fully paired design has both
#' counts.
#'
#' @return A brms custom_family object.
#' @seealso [bipois_partialobs()] for partially observed pairs; [binegbin()]
#'   for the overdispersed counterpart.
#' @export
bipois <- function() {
  brms::custom_family(
    name  = "bipois",
    dpars = c("mu", "lambdaone", "lambdatwo"),  # mu = lambda_shared -- see Details
    links = c("log", "log", "log"),
    lb    = c(0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "1")  # vint1 = y2; the flag is a literal -- see Details
  )
}

#' @rdname bipois
#' @export
bipois_stanvars <- function() {
  brms::stanvar(block = "functions", scode = bipois_stan_funs)
}

#' Joint bivariate-Poisson family for partially observed pairs
#'
#' @description
#' [bipois()] for a design in which the first count is missing on some rows.
#' Same generative model, same `name`, same three dpars, same likelihood and
#' the same post-processing methods -- the only difference is that each row
#' carries a second supplementary integer, an observation flag, through
#' `vint()`:
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
#' full joint [bipois()] lpmf on `(y1, y2)`. A row whose first count was never
#' recorded (`y1_obs == 0`) contributes the second count's marginal *from the
#' same model*. For Poisson components that marginal is closed form -- a sum of
#' independent Poissons is Poisson -- so it is exactly
#' `y2 ~ Poisson(mu + lambdatwo)`. [binegbin_partialobs()] must evaluate the
#' corresponding convolution as a sum; this family does not. Either way the row
#' is not dropped and is not given a different model: it still informs `mu`,
#' `lambdatwo` and any group-level effects.
#'
#' **What you get afterwards.** The fitted model can impute the unobserved
#' first count conditional on the observed second one, which is usually why
#' someone wanted this. `posterior_predict()` and `posterior_epred()` return a
#' `y1` draw and `E[y1 | y2]` for *every* row, matched and unmatched alike --
#' `y1_obs` selects a likelihood branch, not a prediction.
#'
#' **The design consequence, worth knowing before collecting data.**
#' `lambdatwo` appears on both branches, so every row informs it. `lambdaone`
#' appears only on the matched branch and is identified by the matched rows
#' *alone*. `mu` appears on both, but the unmatched branch sees it only through
#' the sum `mu + lambdatwo` -- those rows constrain the total rate of the
#' observed margin, not how it divides between the shared and source-2-only
#' components. Separating `mu` from `lambdatwo`, and so estimating the
#' congruence \eqn{f}, therefore also rests on the matched rows. With few of
#' them, `mu` and `lambdatwo` trade off along their sum and the prior does
#' correspondingly more of the work, however many unmatched rows the design
#' contains.
#'
#' **This is not censoring in brms's sense.** brms's `cens()` addition term
#' means a value known to lie in a set -- `left`, `right`, `interval`. Here the
#' first count is not observed at all and the likelihood marginalises over its
#' whole support. This family was called `bipois_cens()` up to 0.9.1; the name
#' was wrong and was changed at 0.10.0. Do not combine this family with
#' `cens()`.
#'
#' Use in a brm() call as:
#'   brm(
#'     bf(y1 | vint(y2, y1_obs) ~ 1,
#'        mu ~ 1 + (1 | vessel) + (1 | vessel:trip_id),
#'        nlf(lambdaone ~ lamx + methd),
#'        nlf(lambdatwo ~ lamx - methd),
#'        lamx ~ 1, methd ~ 1, nl = TRUE),
#'     family   = bipois_partialobs(),
#'     stanvars = bipois_partialobs_stanvars(),
#'     data     = dat
#'   )
#'
#' @details
#' **When to use this rather than [binegbin_partialobs()].** This family fixes
#' each latent component's variance equal to its mean. Where the counts are
#' genuinely overdispersed relative to that,
#' [binegbin_partialobs()] is the correct model and this one will understate
#' the marginal variances. Where they are not, [binegbin_partialobs()] can only
#' represent the fit by driving its dispersions to their Poisson limit
#' (`shape` \eqn{\to\infty}, equivalently `kappa` \eqn{\to 0}), a boundary at
#' which sampling degrades; fitting the equidispersed family directly avoids
#' it. Compare the two with `loo()`.
#'
#' **One likelihood, two constructors.** This returns the same
#' `custom_family` `name` as [bipois()], so both land on one `bipois_lpmf` and
#' one set of `log_lik_bipois()` / `posterior_predict_bipois()` /
#' `posterior_epred_bipois()` methods. The matched branch is therefore not a
#' second copy that can drift from the fully paired likelihood; it is the same
#' code. See [binegbin_partialobs()] for why `vars` carries a literal in the
#' plain constructor rather than the two families declaring overloaded Stan
#' functions.
#'
#' **Two `vint()` arguments, in declared order.** brms appends `vint()` integers
#' to the generated lpmf call in the order they are listed in the formula's
#' `vint()` term, matching the `vars` declared here
#' (`c("vint1[n]", "vint2[n]")`): so `vint(y2, y1_obs)` binds `vint1 = y2` and
#' `vint2 = y1_obs`. brms generates `target += bipois_lpmf(Y[n] | mu[n],
#' lambdaone[n], lambdatwo[n], vint1[n], vint2[n])`. Reordering the dpars or the
#' two `vint()` terms without matching the Stan signature silently swaps which
#' rate governs which component or which integer is the branch flag.
#'
#' **Telling which shape a stored fit used.** `family$name` is `"bipois"`
#' either way. What distinguishes them is the presence of the second
#' supplementary integer:
#'
#' ```
#' "vint2" %in% names(brms::standata(fit))   # TRUE for a partially observed fit
#' fit$family$vars     # c("vint1[n]", "vint2[n]") or c("vint1[n]", "1")
#' ```
#'
#' @return A brms custom_family object.
#' @seealso [bipois()] for the fully paired case; [binegbin_partialobs()] for
#'   the overdispersed counterpart.
#' @export
bipois_partialobs <- function() {
  brms::custom_family(
    name  = "bipois",
    dpars = c("mu", "lambdaone", "lambdatwo"),  # mu = lambda_shared -- see bipois()
    links = c("log", "log", "log"),
    lb    = c(0, 0, 0),
    type  = "int",
    vars  = c("vint1[n]", "vint2[n]")  # vint1 = y2, vint2 = y1_obs
  )
}

#' @rdname bipois_partialobs
#' @export
bipois_partialobs_stanvars <- function() {
  # Deliberately the same object as bipois_stanvars(): one family name means
  # one Stan function, and both constructors need exactly it. The alias exists
  # so that a call reads consistently, not because there is a second
  # implementation.
  bipois_stanvars()
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
#
# The y1_obs == 0 branch is the closed-form marginal,
# poisson_lpmf(y2 | mu + lambdatwo). No sum, so nothing to cap or truncate;
# this is the branch that has no counterpart in binegbin, where the same
# quantity is a convolution over k = 0..y2.
#
# bipois() reaches this function with y1_obs supplied as the literal 1, so a
# fully paired model always takes the first branch.
bipois_stan_funs <- "
  real bipois_lpmf(int y1, real mu, real lambdaone, real lambdatwo,
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
# common length), used to power log_lik_bipois() and to validate the Stan
# lpmf. Internal reference only, exactly the role binegbin_lpmf_r plays for
# its own family.
#
# The matched branch is computed term-by-term from the original
# P(N_shared=k) P(N1=x-k) P(N2=y-k) definition rather than via the recurrence
# -- an independent route to the same quantity, which is what makes it a check
# on the Stan implementation. It is evaluated post-hoc, not inside the
# sampler's hot loop, so there is no reason to use the recurrence's algebraic
# shortcuts here.
#
# The unmatched branch uses the closed form, as Stan does. The independent
# route -- brute-force convolution over k = 0..y2 -- lives in
# test-bipois-partialobs.R rather than here, so log_lik() does not pay for a
# sum that has an exact one-line answer. Keeping it in the test is what makes
# the closed form checked rather than assumed.
#
# `y1_obs` defaults to 1, so a call that names no flag is the matched branch.
bipois_lpmf_r <- function(y1, y2, mu, lambdaone, lambdatwo, y1_obs = 1L) {
  n <- max(length(y1), length(y2), length(y1_obs), length(mu),
           length(lambdaone), length(lambdatwo))
  y1        <- rep_len(y1, n)
  y2        <- rep_len(y2, n)
  y1_obs    <- rep_len(y1_obs, n)
  mu        <- rep_len(mu, n)
  lambdaone <- rep_len(lambdaone, n)
  lambdatwo <- rep_len(lambdatwo, n)

  vapply(seq_len(n), function(i) {
    if (y1_obs[i] != 1) {
      return(stats::dpois(y2[i], mu[i] + lambdatwo[i], log = TRUE))
    }
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
  y1     <- prep$data$Y[i]
  y2     <- prep$data$vint1[i]
  y1_obs <- .y1_obs_at(prep, i)
  bipois_lpmf_r(y1, y2, mu, lambdaone, lambdatwo, y1_obs)
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
  #
  # y1_obs is deliberately NOT read: every row gets a y1 draw conditional on
  # its observed y2, matched and y2-only alike, which is what imputing the
  # unobserved margin across every row requires.
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
  # expectation rather than simulated. y1_obs is not read, for the same reason
  # it is not read there.
  y2_mat <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}
