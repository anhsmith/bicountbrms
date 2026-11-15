# tests/testthat/test-impute.R
#
# posterior_predict_y1_given_y2() and posterior_epred_y1_given_y2() return a
# draw of y1 given the recorded y2, and E[y1 | y2] = E[N_shared | y2] +
# lambdaone. posterior_predict_y1_and_y2() and posterior_epred_y1_and_y2()
# return draws of the pair (y1, y2) and the marginal means of the two counts.
# The tests check that (1) with keep_observed = FALSE, the conditional
# functions reproduce the reference methods below exactly; (2) keep_observed =
# TRUE returns the recorded y1 on paired rows and leaves unpaired rows
# unchanged; (3) the conditional expectation equals the mean of the
# conditional draws; and (4) the draws of the pair have the marginal moments,
# including Cov(y1, y2) = Var(N_shared). The Stan block at the end fits one
# small model and checks the exported functions through brms.

# ---------------------------------------------------------------------------
# Reference methods for the conditional draw and mean
# ---------------------------------------------------------------------------
#
# These are the bodies of posterior_predict_<family>() and
# posterior_epred_<family>() in 0.10.1, which conditioned on y2.

old_pp_bipois <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  y2 <- prep$data$vint1[i]
  p_shared <- mu / (mu + lambdatwo)
  n_shared <- stats::rbinom(length(mu), size = y2, prob = p_shared)
  n1       <- stats::rpois(length(mu), lambdaone)
  n_shared + n1
}

old_ep_bipois <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  y2 <- prep$data$vint1
  y2_mat <- matrix(y2, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  p_shared <- mu / (mu + lambdatwo)
  y2_mat * p_shared + lambdaone
}

old_pp_binegbin <- function(i, prep, ...) {
  mu        <- brms::get_dpar(prep, "mu", i = i)
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem", i = i)
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb", i = i)
  shapes    <- brms::get_dpar(prep, "shapes", i = i)
  shapexone <- .get_dpar_any(prep, .SHAPEXONE_NAMES, i = i)
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES, i = i)
  y2 <- prep$data$vint1[i]
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

old_ep_binegbin <- function(prep) {
  mu        <- brms::get_dpar(prep, "mu")
  lambdaone <- .get_rate(prep, "lambdaone", "lambdaem")
  lambdatwo <- .get_rate(prep, "lambdatwo", "lambdalb")
  shapes    <- brms::get_dpar(prep, "shapes")
  shapextwo <- .get_dpar_any(prep, .SHAPEXTWO_NAMES)
  y2 <- prep$data$vint1
  .e_shared_given_y2_nb(mu, lambdatwo, shapes, shapextwo, y2) + lambdaone
}

old_pp_matrix <- function(fun, prep) {
  do.call(cbind, lapply(seq_len(prep$nobs), fun, prep = prep))
}

# ---------------------------------------------------------------------------
# A multi-row prep: dpars varying by draw and row, paired and unpaired rows
# ---------------------------------------------------------------------------

NDR <- 300L
Y   <- c(4L, 0L, 0L, 11L, 0L, 7L)
V1  <- c(0L, 3L, 6L, 9L, 2L, 5L)    # y2
V2  <- c(1L, 0L, 0L, 1L, 0L, 1L)    # y1_obs

multi_dpars <- function(nb, seed = 3L) {
  set.seed(seed)
  n <- length(Y)
  m <- function(lo, hi) matrix(stats::runif(NDR * n, lo, hi), NDR, n)
  dp <- list(mu = m(2, 6), lambdaone = m(1, 4), lambdatwo = m(1, 5))
  if (nb) {
    dp$shapes    <- m(0.5, 5)
    dp$shapexone <- m(0.3, 3)
    dp$shapextwo <- m(2, 10)
  }
  dp
}

multi_prep <- function(nb, vint2 = V2) {
  make_synthetic_prep(multi_dpars(nb), Y = Y, vint1 = V1, vint2 = vint2)
}

FAMILIES <- list(
  bipois   = list(nb = FALSE, pp = old_pp_bipois,   ep = old_ep_bipois),
  binegbin = list(nb = TRUE,  pp = old_pp_binegbin, ep = old_ep_binegbin)
)

# ---------------------------------------------------------------------------
# 1. Exact reproduction of the reference methods
# ---------------------------------------------------------------------------

test_that("posterior_predict_y1_given_y2() reproduces the 0.10.1 posterior_predict draws", {
  for (fam in names(FAMILIES)) {
    f <- FAMILIES[[fam]]
    for (v2 in list(V2, NULL)) {
      prep <- multi_prep(f$nb, vint2 = v2)
      set.seed(101); old <- old_pp_matrix(f$pp, prep)
      set.seed(101); new <- .impute_y1_prep(prep, fam, type = "draw",
                                            keep_observed = FALSE)
      expect_identical(new, old, label = paste(fam, "draws"))
    }
  }
})

test_that("posterior_epred_y1_given_y2() reproduces the 0.10.1 posterior_epred", {
  for (fam in names(FAMILIES)) {
    f <- FAMILIES[[fam]]
    prep <- multi_prep(f$nb)
    expect_equal(
      .impute_y1_prep(prep, fam, type = "mean", keep_observed = FALSE),
      f$ep(prep),
      tolerance = 1e-14,
      label = paste(fam, "expectation")
    )
  }
})

# ---------------------------------------------------------------------------
# 2. keep_observed
# ---------------------------------------------------------------------------

test_that("keep_observed = TRUE returns the recorded y1 on paired rows only", {
  paired <- V2 == 1L
  for (fam in names(FAMILIES)) {
    prep <- multi_prep(FAMILIES[[fam]]$nb)
    for (type in c("draw", "mean")) {
      set.seed(5); kept <- .impute_y1_prep(prep, fam, type, keep_observed = TRUE)
      set.seed(5); all  <- .impute_y1_prep(prep, fam, type, keep_observed = FALSE)
      expect_equal(dim(kept), c(NDR, length(Y)))
      expect_equal(kept[, paired],
                   matrix(Y[paired], NDR, sum(paired), byrow = TRUE),
                   label = paste(fam, type, "paired rows"))
      # The unpaired rows are unchanged, draw for draw.
      expect_identical(kept[, !paired], all[, !paired],
                       label = paste(fam, type, "unpaired rows"))
      # On paired rows, the conditional draw or expectation differs from the
      # recorded y1, so the check above is not vacuous.
      expect_false(isTRUE(all.equal(kept[, paired], all[, paired])))
    }
  }
})

test_that("keep_observed = TRUE on a fully paired prep returns the recorded y1 everywhere", {
  # A prep from bipois() or binegbin() has no vint2, and every row is paired.
  for (fam in names(FAMILIES)) {
    prep <- multi_prep(FAMILIES[[fam]]$nb, vint2 = NULL)
    out <- .impute_y1_prep(prep, fam, "draw", keep_observed = TRUE)
    expect_equal(out, matrix(Y, NDR, length(Y), byrow = TRUE))
  }
})

test_that("keep_observed = TRUE stops if a paired row has no recorded y1", {
  # brms fills an omitted response with NA when check_response = FALSE.
  prep <- multi_prep(FALSE)
  prep$data$Y <- rep(NA_integer_, length(Y))
  expect_error(.impute_y1_prep(prep, "bipois", "draw", TRUE),
               "missing on 3 row")
  expect_equal(dim(.impute_y1_prep(prep, "bipois", "draw", FALSE)),
               c(NDR, length(Y)))
})

test_that("a family from outside the package is refused", {
  prep <- multi_prep(FALSE)
  expect_error(.impute_y1_prep(prep, "poisson", "draw", TRUE), "bicountbrms")
  expect_error(.joint_prep(prep, "skellam", "draw"), "bicountbrms")
})

# ---------------------------------------------------------------------------
# 3. The conditional expectation
# ---------------------------------------------------------------------------

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2; SHAPEX <- 1.5
Y2 <- 9L
ND <- 5e4L   # MC standard error on the mean is ~0.02 at these settings

pois_prep <- function(vint2 = NULL, y2 = Y2, ndraws = ND) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ndraws), lambdaone = rep(LONE, ndraws),
                 lambdatwo = rep(LTWO, ndraws)),
    Y = 0L, vint1 = y2, vint2 = vint2
  )
}

nb_prep <- function(vint2 = NULL, y2 = Y2, ndraws = ND, shapex = SHAPEX) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ndraws), lambdaone = rep(LONE, ndraws),
                 lambdatwo = rep(LTWO, ndraws), shapes = rep(SHAPES, ndraws),
                 shapex = rep(shapex, ndraws)),
    Y = 0L, vint1 = y2, vint2 = vint2
  )
}

test_that("the bipois conditional expectation equals the Binomial-split closed form", {
  # Conditioning a sum of independent Poissons on its total gives a Binomial,
  # so E[N_shared | y2] = y2 * mu/(mu + lambdatwo) with no sum to evaluate.
  closed <- Y2 * MU / (MU + LTWO) + LONE
  for (v2 in list(NULL, 1L, 0L)) {
    expect_equal(unique(as.vector(.impute_mean_bipois(pois_prep(v2, ndraws = 2L)))),
                 closed, tolerance = 1e-12)
  }
})

test_that("the binegbin conditional expectation equals the mean of the conditional draws", {
  # The mean of the conditional draws is a reference that assumes no closed
  # form.
  set.seed(20260804)
  draws <- .impute_draw_binegbin(1, nb_prep())
  ep    <- unique(as.vector(.impute_mean_binegbin(nb_prep(ndraws = 2L))))
  expect_length(ep, 1L)
  expect_equal(mean(draws), ep, tolerance = 0.1)

  # The Poisson shared fraction mu/(mu + lambdatwo) must not reproduce the
  # binegbin value, or the check above would be vacuous.
  approx <- Y2 * MU / (MU + LTWO) + LONE
  expect_gt(abs(approx - ep), 0.01)
})

test_that("the binegbin conditional expectation tends to the bipois one as the dispersions grow", {
  # NB2(m, phi) -> Poisson(m) as phi -> Inf, for shapes and shapex together.
  ref <- .impute_mean_bipois(pois_prep(ndraws = 2L))
  err <- function(phi) {
    p <- make_synthetic_prep(
      dpars = list(mu = rep(MU, 2L), lambdaone = rep(LONE, 2L),
                   lambdatwo = rep(LTWO, 2L), shapes = rep(phi, 2L),
                   shapex = rep(phi, 2L)),
      Y = 0L, vint1 = Y2
    )
    max(abs(.impute_mean_binegbin(p) - ref))
  }
  e5 <- err(1e5)
  e7 <- err(1e7)
  expect_lt(e5, 1e-2)
  expect_lt(e7, 1e-4)
  expect_lt(e7, e5)          # the error decreases as phi increases
})

test_that("y2 = 0 leaves the conditional expectation at lambdaone alone", {
  # N_shared | y2 = 0 is degenerate at 0, so E[y1 | y2 = 0] = lambdaone exactly.
  expect_equal(unique(as.vector(.impute_mean_bipois(pois_prep(1L, y2 = 0L, ndraws = 2L)))),
               LONE, tolerance = 1e-12)
  expect_equal(unique(as.vector(.impute_mean_binegbin(nb_prep(1L, y2 = 0L, ndraws = 2L)))),
               LONE, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. posterior_predict_y1_and_y2(): marginal moments
# ---------------------------------------------------------------------------

SX1 <- 0.6; SX2 <- 7
NP  <- 2e5L

pair_prep <- function(nb) {
  dp <- list(mu = rep(MU, NP), lambdaone = rep(LONE, NP),
             lambdatwo = rep(LTWO, NP))
  if (nb) {
    dp$shapes    <- rep(SHAPES, NP)
    dp$shapexone <- rep(SX1, NP)
    dp$shapextwo <- rep(SX2, NP)
  }
  make_synthetic_prep(dp, Y = 0L, vint1 = Y2, vint2 = 0L)
}

nb_var <- function(m, phi) m + m^2 / phi

test_that("the bipois pair draws have the analytic marginal moments", {
  set.seed(20261114)
  out <- .joint_prep(pair_prep(FALSE), "bipois", "draw")
  expect_equal(dim(out), c(NP, 1L, 2L))
  expect_identical(dimnames(out)[[3]], c("y1", "y2"))
  y1 <- out[, 1, "y1"]; y2 <- out[, 1, "y2"]
  expect_equal(mean(y1), MU + LONE, tolerance = 0.01)
  expect_equal(mean(y2), MU + LTWO, tolerance = 0.01)
  expect_equal(var(y1),  MU + LONE, tolerance = 0.02)
  expect_equal(var(y2),  MU + LTWO, tolerance = 0.02)
  expect_equal(cov(y1, y2), MU, tolerance = 0.03)
})

test_that("the binegbin pair draws have the analytic marginal moments", {
  set.seed(20261114)
  out <- .joint_prep(pair_prep(TRUE), "binegbin", "draw")
  y1 <- out[, 1, "y1"]; y2 <- out[, 1, "y2"]
  var_s <- nb_var(MU, SHAPES)
  expect_equal(mean(y1), MU + LONE, tolerance = 0.01)
  expect_equal(mean(y2), MU + LTWO, tolerance = 0.01)
  expect_equal(var(y1), var_s + nb_var(LONE, SX1), tolerance = 0.05)
  expect_equal(var(y2), var_s + nb_var(LTWO, SX2), tolerance = 0.05)
  expect_equal(cov(y1, y2), var_s, tolerance = 0.05)
  # Var(d) = Var(N1) + Var(N2), because the shared component cancels.
  expect_equal(var(y1 - y2), nb_var(LONE, SX1) + nb_var(LTWO, SX2),
               tolerance = 0.05)
  # Exchanging the two exclusive dispersions would change Var(y1) from 35.5 to
  # 21.8, far outside the 5% tolerance.
  expect_gt(abs(nb_var(LONE, SX2) - nb_var(LONE, SX1)) / (var_s + nb_var(LONE, SX1)), 0.3)
})

test_that("the pair means are mu + lambdaone and mu + lambdatwo, and the first equals the epred", {
  for (fam in names(FAMILIES)) {
    prep <- multi_prep(FAMILIES[[fam]]$nb)
    out  <- .joint_prep(prep, fam, "mean")
    expect_equal(dim(out), c(NDR, length(Y), 2L))
    expect_equal(out[, , "y1"], prep$dpars$mu + prep$dpars$lambdaone)
    expect_equal(out[, , "y2"], prep$dpars$mu + prep$dpars$lambdatwo)
    epred <- get(paste0("posterior_epred_", fam))(prep)
    expect_equal(out[, , "y1"], epred)
  }
})

test_that("the pair draws do not read the recorded counts or the indicator", {
  for (fam in names(FAMILIES)) {
    a_prep <- multi_prep(FAMILIES[[fam]]$nb)
    b_prep <- multi_prep(FAMILIES[[fam]]$nb, vint2 = NULL)
    b_prep$data$vint1 <- rev(V1)
    b_prep$data$Y     <- rev(Y)
    set.seed(9); a <- .joint_prep(a_prep, fam, "draw")
    set.seed(9); b <- .joint_prep(b_prep, fam, "draw")
    expect_identical(a, b, label = fam)
  }
})

# ---------------------------------------------------------------------------
# Against a fitted model
# ---------------------------------------------------------------------------

test_that("the _y1_given_y2 and _y1_and_y2 functions work on a fit, with and without newdata", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_no_stan()

  set.seed(20261114)
  n   <- 60L
  x   <- stats::rnorm(n)
  ns  <- stats::rpois(n, exp(1 + 0.4 * x))
  y1  <- ns + stats::rpois(n, 1.5)
  y2  <- ns + stats::rpois(n, 2)
  obs <- rep(c(1L, 0L, 0L), length.out = n)
  dat <- data.frame(x, y1 = ifelse(obs == 1L, y1, 0L), y2, y1_obs = obs)

  fit <- suppressWarnings(suppressMessages(brms::brm(
    brms::bf(y1 | vint(y2, y1_obs) ~ x),
    family   = bipois_partialobs(),
    stanvars = bipois_partialobs_stanvars(),
    data     = dat,
    chains   = 1,
    iter     = 400,
    refresh  = 0,
    seed     = 1
  )))

  # A copy of the fit holds the reference methods in its family object, and
  # brms calls these in preference to the methods of the package
  # (brms:::custom_family_method()).
  old <- fit
  old$family$posterior_predict   <- old_pp_bipois
  old$family$posterior_epred     <- old_ep_bipois
  old$formula$family             <- old$family
  # Both calls of each comparison receive the same draw_ids, because ndraws
  # alone selects a random subset of posterior draws on each call.
  ids <- seq_len(100)
  set.seed(77); pp_old <- brms::posterior_predict(old, draw_ids = ids)
  set.seed(77); pp_new <- posterior_predict_y1_given_y2(fit, draw_ids = ids)
  # brms attaches an empty dimnames list, so the values are compared unnamed.
  expect_equal(unname(pp_new), unname(pp_old))
  ep_old <- brms::posterior_epred(old, draw_ids = ids)
  ep_new <- posterior_epred_y1_given_y2(fit, draw_ids = ids)
  expect_equal(unname(ep_new), unname(ep_old))

  # posterior_epred() on the original fit returns the marginal mean, which
  # differs from the reference, so the copy did dispatch to the reference
  # methods.
  ep_marg <- brms::posterior_epred(fit, draw_ids = ids)
  expect_false(isTRUE(all.equal(unname(ep_marg), unname(ep_old))))
  prep <- brms::prepare_predictions(fit, draw_ids = ids)
  expect_equal(unname(ep_marg),
               unname(brms::get_dpar(prep, "mu") +
                        brms::get_dpar(prep, "lambdaone")))

  # keep_observed = TRUE on the data of the fit.
  imp <- posterior_predict_y1_given_y2(fit, keep_observed = TRUE, ndraws = 50)
  expect_equal(imp[, obs == 1L],
               matrix(dat$y1[obs == 1L], 50, sum(obs), byrow = TRUE))

  # newdata with real y2: only the rows supplied.
  nd  <- dat[1:5, ]
  out <- posterior_predict_y1_given_y2(fit, newdata = nd, ndraws = 20)
  expect_equal(dim(out), c(20L, 5L))

  # newdata without y1: brms requires the vint() variables but not the
  # response. keep_observed = TRUE then requires y1 on the paired rows.
  no_y1 <- nd[, c("x", "y2", "y1_obs")]
  expect_equal(
    dim(posterior_predict_y1_given_y2(fit, newdata = no_y1, ndraws = 20)),
    c(20L, 5L)
  )
  expect_error(
    posterior_predict_y1_given_y2(fit, newdata = no_y1, ndraws = 20,
                                    keep_observed = TRUE),
    "keep_observed"
  )
  unpaired_only <- transform(no_y1, y1_obs = 0L)
  expect_equal(
    dim(posterior_predict_y1_given_y2(fit, newdata = unpaired_only,
                                        ndraws = 20, keep_observed = TRUE)),
    c(20L, 5L)
  )
  # The default is keep_observed = FALSE.
  expect_identical(formals(posterior_predict_y1_given_y2)$keep_observed, FALSE)
  expect_identical(formals(posterior_epred_y1_given_y2)$keep_observed, FALSE)
  expect_equal(dim(posterior_predict_y1_and_y2(fit, newdata = no_y1, ndraws = 20)),
               c(20L, 5L, 2L))
  # Without y2 or the observation indicator, brms stops with an error.
  expect_error(posterior_predict_y1_and_y2(fit, newdata = nd[, "x", drop = FALSE]))

  pairs <- posterior_predict_y1_and_y2(fit, ndraws = 30)
  expect_equal(dim(pairs), c(30L, n, 2L))
  expect_true(all(pairs >= 0))
  means <- posterior_epred_y1_and_y2(fit, draw_ids = 1:30)
  expect_equal(unname(means[, , "y1"]),
               unname(brms::posterior_epred(fit, draw_ids = 1:30)))

  # newdata with placeholder counts: the marginal quantities do not read them.
  ph <- transform(nd, y1 = 0L, y2 = 0L, y1_obs = 0L)
  expect_equal(posterior_epred_y1_and_y2(fit, newdata = ph, draw_ids = 1:10),
               posterior_epred_y1_and_y2(fit, newdata = nd, draw_ids = 1:10))
  expect_equal(brms::posterior_epred(fit, newdata = ph, draw_ids = 1:10),
               brms::posterior_epred(fit, newdata = nd, draw_ids = 1:10))
})
