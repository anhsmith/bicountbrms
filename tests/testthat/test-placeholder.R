# tests/testthat/test-placeholder.R
#
# The placeholder for y1 on unpaired rows. brms drops rows with a missing
# response, so a fit from bipois_partialobs() or binegbin_partialobs()
# requires a value of y1 on every row. On an unpaired row, the likelihood term
# is the marginal probability of y2, which does not depend on y1. These tests
# check that the R references and every post-processing method give identical
# results for the placeholders 0, -1, and -5. No Stan compilation is required.

PH <- c(0L, -1L, -5L)
y2_grid <- c(0L, 1L, 4L, 12L)

test_that("bipois_lpmf_r() on unpaired rows does not depend on the placeholder", {
  for (y2 in y2_grid) {
    vals <- vapply(PH, function(ph) {
      bipois_lpmf_r(ph, y2, mu = 3.2, lambdaone = 1.4, lambdatwo = 2.1,
                    y1_obs = 0L)
    }, numeric(1))
    expect_identical(vals[2], vals[1], label = paste("y2 =", y2, "placeholder -1"))
    expect_identical(vals[3], vals[1], label = paste("y2 =", y2, "placeholder -5"))
    expect_true(is.finite(vals[1]))
  }
})

test_that("binegbin_lpmf_r() on unpaired rows does not depend on the placeholder", {
  for (y2 in y2_grid) {
    vals <- vapply(PH, function(ph) {
      binegbin_lpmf_r(ph, y2, mu = 3.2, lambdaone = 1.4, lambdatwo = 2.1,
                      shapes = 2, shapexone = 0.7, shapextwo = 5,
                      y1_obs = 0L)
    }, numeric(1))
    expect_identical(vals[2], vals[1], label = paste("y2 =", y2, "placeholder -1"))
    expect_identical(vals[3], vals[1], label = paste("y2 =", y2, "placeholder -5"))
    expect_true(is.finite(vals[1]))
  }
})

# A prep with three paired and three unpaired rows, built with the placeholder
# `ph` on the unpaired rows.
ph_prep <- function(nb, ph) {
  y1  <- c(4L, 0L, 7L, 0L, 2L, 0L)
  obs <- c(1L, 0L, 1L, 0L, 1L, 0L)
  y1[obs == 0L] <- ph
  dp <- list(mu = c(3, 3.5, 2.5), lambdaone = c(1.2, 1.5, 1),
             lambdatwo = c(2, 2.4, 1.8))
  if (nb) {
    dp$shapes <- c(2, 3, 1.5)
    dp$shapexone <- c(0.8, 1, 0.6)
    dp$shapextwo <- c(4, 5, 3)
  }
  make_synthetic_prep(dp, Y = y1, vint1 = c(5L, 3L, 9L, 0L, 4L, 6L),
                      vint2 = obs)
}

test_that("log_lik and the prediction functions do not depend on the placeholder", {
  for (fam in c("bipois", "binegbin")) {
    nb <- fam == "binegbin"
    z <- ph_prep(nb, 0L)
    m <- ph_prep(nb, -1L)
    ll <- get(paste0("log_lik_", fam))
    pp <- get(paste0("posterior_predict_", fam))
    ep <- get(paste0("posterior_epred_", fam))
    for (i in seq_len(z$nobs)) {
      expect_identical(ll(i, m), ll(i, z), label = paste(fam, "log_lik, obs", i))
      set.seed(3); a <- pp(i, m)
      set.seed(3); b <- pp(i, z)
      expect_identical(a, b, label = paste(fam, "posterior_predict, obs", i))
    }
    expect_identical(ep(m), ep(z))
    for (type in c("draw", "mean")) {
      for (keep in c(TRUE, FALSE)) {
        set.seed(3); a <- .impute_y1_prep(m, fam, type, keep_observed = keep)
        set.seed(3); b <- .impute_y1_prep(z, fam, type, keep_observed = keep)
        expect_identical(a, b, label = paste(fam, "imputation", type, keep))
      }
      set.seed(3); a <- .joint_prep(m, fam, type)
      set.seed(3); b <- .joint_prep(z, fam, type)
      expect_identical(a, b, label = paste(fam, "joint", type))
    }
  }
})
