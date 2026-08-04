# tests/testthat/test-deprecated.R
#
# The 0.9.0 rename binegbin_joint -> binegbin_cens leaves five names behind.
# Two of them (the constructor and stanvars) are a convenience; the other three
# are load-bearing, because a fit made before 0.9.0 stores
# family$name == "binegbin_joint" forever and brms builds log_lik_<name> /
# posterior_predict_<name> / posterior_epred_<name> from that stored name. If
# any of the three disappears, loo() and posterior_predict() on a stored fit
# fail with "could not find function" and nothing the caller passes can route
# around it.
#
# These tests therefore assert two things that are easy to break silently:
# that all five names are exported, and that the three post-processing shims
# are value-identical to their new counterparts. R-side only.

MU <- 5; LONE <- 3; LTWO <- 4; SHAPES <- 2; SHAPEX1 <- 1.5; SHAPEX2 <- 2.5
Y1 <- 6L; Y2 <- 9L

dep_prep <- function(ndraws = 200L, y1_obs = 1L) {
  make_synthetic_prep(
    dpars = list(mu = rep(MU, ndraws), lambdaone = rep(LONE, ndraws),
                 lambdatwo = rep(LTWO, ndraws), shapes = rep(SHAPES, ndraws),
                 shapexone = rep(SHAPEX1, ndraws), shapextwo = rep(SHAPEX2, ndraws)),
    Y = Y1, vint1 = Y2, vint2 = y1_obs
  )
}

test_that("every deprecated name is still exported", {
  # Exported, not merely present: brms resolves these off the search path, so
  # an internal-only function would not be found from a user's session.
  exported <- getNamespaceExports("bicountbrms")
  for (nm in c("binegbin_joint", "binegbin_joint_stanvars",
               "log_lik_binegbin_joint", "posterior_predict_binegbin_joint",
               "posterior_epred_binegbin_joint")) {
    expect_true(nm %in% exported, label = nm)
  }
})

test_that("the deprecated constructor warns and returns the new family", {
  expect_warning(fam <- binegbin_joint(), "deprecated")
  # It returns binegbin_cens itself, not a copy carrying the old name. A fit
  # made through the deprecated constructor is therefore a binegbin_cens fit
  # and needs no shims of its own -- the compatibility problem does not
  # propagate forward.
  expect_identical(fam, binegbin_cens())
  expect_identical(fam$name, "binegbin_cens")

  expect_warning(sv <- binegbin_joint_stanvars(), "deprecated")
  expect_identical(sv, binegbin_cens_stanvars())
})

test_that("the post-processing shims forward without warning", {
  # Silence is the specification, not an oversight: brms calls the first two
  # once per observation, and the stored family name is not something the
  # caller can change.
  prep <- dep_prep()
  expect_silent(ll <- log_lik_binegbin_joint(1, prep))
  expect_equal(ll, log_lik_binegbin_cens(1, prep))

  expect_silent(ep <- posterior_epred_binegbin_joint(prep))
  expect_equal(ep, posterior_epred_binegbin_cens(prep))

  set.seed(11)
  pp_old <- posterior_predict_binegbin_joint(1, dep_prep())
  set.seed(11)
  pp_new <- posterior_predict_binegbin_cens(1, dep_prep())
  expect_identical(pp_old, pp_new)
})

test_that("the shims forward on the censored branch too", {
  # The branch that only the censoring-aware family has, and the one a stored
  # fit is most likely to exercise.
  prep <- dep_prep(y1_obs = 0L)
  expect_equal(log_lik_binegbin_joint(1, prep), log_lik_binegbin_cens(1, prep))
  expect_equal(posterior_epred_binegbin_joint(prep),
               posterior_epred_binegbin_cens(prep))
})

test_that("pre-0.8.0 five-dpar fits still resolve through the deprecated names", {
  # The two compatibility layers compose: a fit made at 0.7.0 spells its single
  # excess dispersion `shapex` AND names its family `binegbin_joint`. It has to
  # pass through the shim here and .get_dpar_any()'s `shapex` fallback in
  # utils.R, and the result must equal the six-dpar answer with the two excess
  # dispersions tied.
  five <- make_synthetic_prep(
    dpars = list(mu = rep(MU, 50L), lambdaone = rep(LONE, 50L),
                 lambdatwo = rep(LTWO, 50L), shapes = rep(SHAPES, 50L),
                 shapex = rep(SHAPEX1, 50L)),
    Y = Y1, vint1 = Y2, vint2 = 1L
  )
  six <- make_synthetic_prep(
    dpars = list(mu = rep(MU, 50L), lambdaone = rep(LONE, 50L),
                 lambdatwo = rep(LTWO, 50L), shapes = rep(SHAPES, 50L),
                 shapexone = rep(SHAPEX1, 50L), shapextwo = rep(SHAPEX1, 50L)),
    Y = Y1, vint1 = Y2, vint2 = 1L
  )
  expect_equal(log_lik_binegbin_joint(1, five), log_lik_binegbin_cens(1, six))
})

test_that("bipois_cens has no deprecated aliases", {
  # It was added in 0.9.0 under its final name, so no fit and no source can
  # refer to it as bipois_joint(). Shims for it would be dead code that the
  # next major version would have to justify removing.
  exported <- getNamespaceExports("bicountbrms")
  expect_false(any(grepl("bipois_joint", exported)))
})
