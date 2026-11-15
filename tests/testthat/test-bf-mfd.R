# tests/testthat/test-bf-mfd.R
#
# bf_mfd() builds the brms non-linear formula in the congruence-bias
# coordinates (M, f, delta, and kappa). These tests check (i) that the Stan
# code generated from bf_mfd() is identical to the Stan code from the
# hand-written formula, for each of the three valid sets of coordinates, (ii)
# that the rate expressions reproduce mfd_to_dpars(), and (iii) the validation
# of the arguments. stancode() does not compile, so the file requires brms but
# no Stan toolchain.

skip_if_not_installed("brms")

set.seed(20261115)
n  <- 30L
dd <- data.frame(
  y1     = rpois(n, 4),
  y2     = rpois(n, 4),
  y1_obs = rep(c(1L, 0L), length.out = n),
  site   = factor(rep(c("A", "B", "C"), length.out = n))
)

# The hand-written response and rate lines, common to all three sets.
hand_rates <- function() {
  brms::bf(y1 | vint(y2) ~ 1, nl = TRUE) +
    brms::nlf(mu        ~ logM + log_inv_logit(logitf)) +
    brms::nlf(lambdaone ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(2 * delta)) +
    brms::nlf(lambdatwo ~ log(2) + logM + log1m_inv_logit(logitf) + log_inv_logit(-2 * delta))
}

stancode_chr <- function(formula, family, stanvars) {
  as.character(
    brms::stancode(formula, data = dd, family = family, stanvars = stanvars)
  )
}

# ---------------------------------------------------------------------------
# Generated Stan code equals the hand-written formula
# ---------------------------------------------------------------------------

test_that("separate overdispersions give the hand-written binegbin model", {
  hand <- hand_rates() +
    brms::nlf(shapes    ~ -2 * log(kappas)) +
    brms::nlf(shapexone ~ -2 * log(kappaxone)) +
    brms::nlf(shapextwo ~ -2 * log(kappaxtwo)) +
    brms::lf(logM ~ 0 + site, logitf ~ 1, delta ~ 1,
             kappas ~ 1, kappaxone ~ 1, kappaxtwo ~ 1)
  helper <- bf_mfd(
    y1 | vint(y2),
    logM      ~ 0 + site,
    logitf    ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  )
  expect_identical(
    stancode_chr(helper, binegbin(), binegbin_stanvars()),
    stancode_chr(hand,   binegbin(), binegbin_stanvars())
  )
})

test_that("one kappax gives the hand-written tied binegbin model", {
  hand <- hand_rates() +
    brms::nlf(shapes    ~ -2 * log(kappas)) +
    brms::nlf(shapexone ~ -2 * log(kappax)) +
    brms::nlf(shapextwo ~ -2 * log(kappax)) +
    brms::lf(logM ~ 0 + site, logitf ~ 1, delta ~ 1, kappas ~ 1, kappax ~ 1)
  helper <- bf_mfd(
    y1 | vint(y2),
    logM   ~ 0 + site,
    logitf ~ 1,
    delta  ~ 1,
    kappas ~ 1,
    kappax ~ 1
  )
  expect_identical(
    stancode_chr(helper, binegbin(), binegbin_stanvars()),
    stancode_chr(hand,   binegbin(), binegbin_stanvars())
  )
})

test_that("no overdispersion gives the hand-written bipois model", {
  hand <- hand_rates() +
    brms::lf(logM ~ 0 + site, logitf ~ 1, delta ~ 1)
  helper <- bf_mfd(
    y1 | vint(y2),
    logM   ~ 0 + site,
    logitf ~ 1,
    delta  ~ 1
  )
  expect_identical(
    stancode_chr(helper, bipois(), bipois_stanvars()),
    stancode_chr(hand,   bipois(), bipois_stanvars())
  )
})

test_that("the _partialobs response vint(y2, y1_obs) is passed through", {
  helper <- bf_mfd(
    y1 | vint(y2, y1_obs),
    logM   ~ 0 + site,
    logitf ~ 1,
    delta  ~ 1,
    kappas ~ 1,
    kappax ~ 1
  )
  expect_identical(
    deparse(helper$formula),
    "y1 | vint(y2, y1_obs) ~ 1"
  )
  sc <- stancode_chr(
    helper,
    binegbin_partialobs(),
    binegbin_partialobs_stanvars()
  )
  expect_match(sc, "vint1[n], vint2[n])", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# The rate expressions reproduce mfd_to_dpars()
# ---------------------------------------------------------------------------

test_that("the nlf() expressions give the rates and dispersions of mfd_to_dpars()", {
  grid <- expand.grid(
    M     = c(0.5, 3, 12, 40),
    f     = c(0.05, 0.3, 0.67, 0.95),
    delta = c(-1.5, -0.2, 0, 0.4, 2)
  )
  kappas    <- 0.4
  kappaxone <- 0.9
  kappaxtwo <- 0.3

  bform <- bf_mfd(
    y1 | vint(y2),
    logM      ~ 1,
    logitf    ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  )
  # log_inv_logit() and log1m_inv_logit() are the R versions in the
  # brms namespace, which brms uses when it evaluates a non-linear formula
  # in post-processing.
  vals <- list2env(
    list(
      logM      = log(grid$M),
      logitf    = stats::qlogis(grid$f),
      delta     = grid$delta,
      kappas    = kappas,
      kappaxone = kappaxone,
      kappaxtwo = kappaxtwo
    ),
    parent = asNamespace("brms")
  )
  eval_dpar <- function(dpar) {
    exp(eval(bform$pforms[[dpar]][[3]], vals))
  }

  ref <- mfd_to_dpars(
    M         = grid$M,
    f         = grid$f,
    delta     = grid$delta,
    kappas    = kappas,
    kappaxone = kappaxone,
    kappaxtwo = kappaxtwo
  )
  for (dpar in c("mu", "lambdaone", "lambdatwo")) {
    expect_equal(eval_dpar(dpar), ref[[dpar]], tolerance = 1e-12)
  }
  for (dpar in c("shapes", "shapexone", "shapextwo")) {
    expect_equal(
      rep_len(eval_dpar(dpar), nrow(grid)),
      ref[[dpar]],
      tolerance = 1e-12
    )
  }
})

# At delta = 25, tanh(delta) rounds to 1 in double precision, so a rate line
# written with log1m(tanh(delta)) would give -Inf for log(lambdatwo).
test_that("the log rate lines stay exact at an extreme source bias", {
  bform <- bf_mfd(
    y1 | vint(y2),
    logM   ~ 1,
    logitf ~ 1,
    delta  ~ 1
  )
  vals <- list2env(
    list(logM = log(12), logitf = stats::qlogis(0.67), delta = c(-25, 25)),
    parent = asNamespace("brms")
  )
  log_rate <- function(dpar) eval(bform$pforms[[dpar]][[3]], vals)
  base <- log(2) + log(12) + log(1 - 0.67)

  expect_equal(
    log_rate("lambdaone"),
    base + stats::plogis(c(-50, 50), log.p = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    log_rate("lambdatwo"),
    base + stats::plogis(c(50, -50), log.p = TRUE),
    tolerance = 1e-12
  )
})

# ---------------------------------------------------------------------------
# Argument order and validation
# ---------------------------------------------------------------------------

test_that("the order of the coordinate formulae does not matter", {
  canonical <- bf_mfd(
    y1 | vint(y2),
    logM      ~ 0 + site,
    logitf    ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    kappaxone ~ 1,
    kappaxtwo ~ 1
  )
  shuffled <- bf_mfd(
    y1 | vint(y2),
    kappaxtwo ~ 1,
    delta     ~ 1,
    kappas    ~ 1,
    logM      ~ 0 + site,
    kappaxone ~ 1,
    logitf    ~ 1
  )
  expect_identical(names(shuffled$pforms), names(canonical$pforms))
  expect_identical(
    stancode_chr(shuffled,  binegbin(), binegbin_stanvars()),
    stancode_chr(canonical, binegbin(), binegbin_stanvars())
  )
})

test_that("an invalid set of coordinates stops with an error naming the valid sets", {
  # kappax together with kappaxone and kappaxtwo
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logitf ~ 1, delta ~ 1,
           kappas ~ 1, kappax ~ 1, kappaxone ~ 1, kappaxtwo ~ 1),
    "three valid sets"
  )
  # kappas without the exclusive components
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logitf ~ 1, delta ~ 1, kappas ~ 1),
    "logM, logitf, delta, kappas, kappax "
  )
  # delta omitted
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logitf ~ 1),
    "logM, logitf, delta \\(bipois\\)"
  )
})

test_that("malformed arguments stop with an error", {
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logitf ~ 1, delta ~ 1, lambdaone ~ 1),
    "Unknown coordinate `lambdaone`"
  )
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logM ~ site, logitf ~ 1, delta ~ 1),
    "More than one formula for `logM`"
  )
  expect_error(
    bf_mfd(y1 | vint(y2), logM ~ 1, logitf ~ 1, ~ 1),
    "two-sided formula"
  )
  expect_error(
    bf_mfd(logM ~ 1, y1 | vint(y2), logitf ~ 1, delta ~ 1),
    "first argument of `bf_mfd\\(\\)` is the response"
  )
  expect_error(bf_mfd(), "Supply the response")
})
