# tests/testthat/test-stancode-shape.R
#
# THE GENERATED STAN CALL, AND THE LITERAL IN `vars`.
#
# 0.10.0 gives each component distribution one Stan function taking the
# observation flag, and two constructors that differ only in how the flag is
# supplied:
#
#   binegbin()             vars = c("vint1[n]", "1")
#   binegbin_partialobs()  vars = c("vint1[n]", "vint2[n]")
#
# The second entry of the plain constructor's `vars` is not a variable. It is a
# literal that brms pastes into the generated call, so the fully paired model
# reaches the same lpmf with y1_obs fixed at 1 and never asks the user for a
# flag column. The literal is what gives one Stan function instead of two overloaded
# ones, and with it a package that needs no floor on the Stan version -- user
# defined function overloading arrived in Stan 2.29, and DESCRIPTION sets no
# floor on rstan or cmdstanr.
#
# WHY THIS FILE EXISTS. The mechanism depends on brms behaviour that is not a
# documented guarantee. custom_family() validates `vars` no further than
# as.character(), and brms:::stan_log_lik_custom() is its only consumer: it
# strips any "[...]" index, decides whether the entry is an addition term
# needing a response suffix, and pastes the result into the call. A literal
# survives all three steps untouched. If a future brms starts validating
# `vars`, or resolves entries against the data, the failure would otherwise
# surface as an opaque Stan compile error in a user's model rather than here.
#
# This needs brms but no Stan toolchain -- stancode() and standata() do not
# compile -- so it runs in the fast suite.

skip_if_not_installed("brms")

lpmf_call <- function(fam, formula, stanvars, data) {
  sc <- brms::stancode(formula, data = data, family = fam, stanvars = stanvars)
  lines <- strsplit(as.character(sc), "\n")[[1]]
  hit <- grep("target += ", lines, value = TRUE, fixed = TRUE)
  hit <- grep(paste0(fam$name, "_lpmf"), hit, value = TRUE, fixed = TRUE)
  testthat::expect_length(hit, 1L)
  trimws(hit)
}

set.seed(20260825)
n  <- 30L
sd <- data.frame(
  y1     = rpois(n, 4),
  y2     = rpois(n, 4),
  y1_obs = rep(c(1L, 0L), length.out = n)
)

# ---------------------------------------------------------------------------
# The plain constructors pass the flag as a literal
# ---------------------------------------------------------------------------

test_that("binegbin() generates a call ending in the literal 1", {
  call <- lpmf_call(
    binegbin(),
    brms::bf(y1 | vint(y2) ~ 1, lambdaone ~ 1, lambdatwo ~ 1,
             shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
    binegbin_stanvars(), sd
  )
  expect_match(call, "vint1[n], 1)", fixed = TRUE)
  expect_false(grepl("vint2", call))
})

test_that("bipois() generates a call ending in the literal 1", {
  call <- lpmf_call(
    bipois(),
    brms::bf(y1 | vint(y2) ~ 1, lambdaone ~ 1, lambdatwo ~ 1),
    bipois_stanvars(), sd
  )
  expect_match(call, "vint1[n], 1)", fixed = TRUE)
  expect_false(grepl("vint2", call))
})

test_that("the one-vint standata contains no vint2", {
  # The other half of the same fact: the literal exists precisely because
  # there is no second integer column to index.
  sdta <- brms::standata(
    brms::bf(y1 | vint(y2) ~ 1, lambdaone ~ 1, lambdatwo ~ 1,
             shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
    data = sd, family = binegbin(), stanvars = binegbin_stanvars()
  )
  expect_true("vint1" %in% names(sdta))
  expect_false("vint2" %in% names(sdta))
})

# ---------------------------------------------------------------------------
# The _partialobs constructors pass it as data
# ---------------------------------------------------------------------------

test_that("binegbin_partialobs() generates a call ending in vint2[n]", {
  call <- lpmf_call(
    binegbin_partialobs(),
    brms::bf(y1 | vint(y2, y1_obs) ~ 1, lambdaone ~ 1, lambdatwo ~ 1,
             shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
    binegbin_partialobs_stanvars(), sd
  )
  expect_match(call, "vint1[n], vint2[n])", fixed = TRUE)
})

test_that("bipois_partialobs() generates a call ending in vint2[n]", {
  call <- lpmf_call(
    bipois_partialobs(),
    brms::bf(y1 | vint(y2, y1_obs) ~ 1, lambdaone ~ 1, lambdatwo ~ 1),
    bipois_partialobs_stanvars(), sd
  )
  expect_match(call, "vint1[n], vint2[n])", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# One family name, one lpmf, one set of dpars
# ---------------------------------------------------------------------------

test_that("each constructor pair shares a name, dpars and links", {
  # Everything except `vars` must match, because everything except `vars` is
  # what brms uses to build the Stan function name and the parameter block.
  for (pair in list(list(binegbin(), binegbin_partialobs()),
                    list(bipois(),   bipois_partialobs()))) {
    plain <- pair[[1]]; partial <- pair[[2]]
    expect_identical(plain$name,  partial$name)
    expect_identical(plain$dpars, partial$dpars)
    expect_identical(plain$link,  partial$link)
    expect_identical(plain$lb,    partial$lb)
    expect_false(identical(plain$vars, partial$vars))
  }
})

test_that("the stanvars aliases inject exactly the same Stan code", {
  # One likelihood implementation per component distribution is the aim of
  # the unification; the alias is a naming convenience, not a second copy.
  expect_identical(binegbin_stanvars(), binegbin_partialobs_stanvars())
  expect_identical(bipois_stanvars(),   bipois_partialobs_stanvars())
})

test_that("both negative-binomial constructors declare the full six dpars", {
  want <- c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo")
  expect_identical(binegbin()$dpars, want)
  expect_identical(binegbin_partialobs()$dpars, want)
  # And `shapex` is gone: nothing shipping accepts it any more.
  expect_false("shapex" %in% binegbin()$dpars)
})

# ---------------------------------------------------------------------------
# The symmetric special case: tying the two excess dispersions
# ---------------------------------------------------------------------------
#
# binegbin() declares shapexone and shapextwo separately from 0.10.0. The way
# back to the single-dispersion model is a FORMULA constraint, routing both
# through one non-linear parameter:
#
#   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1
#
# That recipe is documented in ?binegbin, in NEWS, in the get-started vignette
# and in two articles, so it is worth a guard. Two things can silently go
# wrong and neither shows up as an error:
#
#   * the tie not actually tying -- if brms generated two independent
#     parameters the model would still compile, still sample, and quietly fit
#     the six-dpar model the user was trying to constrain;
#   * the prior spelling -- routing a dpar through a non-linear parameter moves
#     its prior from class = "Intercept", dpar = "shapex" to class = "b",
#     nlpar = "shapexx". BOTH fields change. A prior written the old way is
#     silently dropped, leaving the parameter flat and improper.
#
# Needs brms but no Stan toolchain, so this runs in the fast suite.

tied_formula <- function() {
  brms::bf(y1 | vint(y2) ~ 1, mu ~ 1,
           brms::nlf(lambdaone ~ lamx), brms::nlf(lambdatwo ~ lamx),
           brms::nlf(shapexone ~ shapexx), brms::nlf(shapextwo ~ shapexx),
           lamx ~ 1, shapes ~ 1, shapexx ~ 1, nl = TRUE)
}

test_that("nlf tying really gives both dispersions one parameter", {
  sc <- brms::stancode(tied_formula(), data = sd, family = binegbin(),
                       stanvars = binegbin_stanvars())
  lines <- strsplit(as.character(sc), "\n")[[1]]

  one <- trimws(grep("shapexone[n] =", lines, fixed = TRUE, value = TRUE))
  two <- trimws(grep("shapextwo[n] =", lines, fixed = TRUE, value = TRUE))
  expect_length(one, 1L)
  expect_length(two, 1L)

  # The substantive check: both are assigned from the SAME non-linear
  # parameter, so there is one free quantity rather than two.
  expect_match(one, "nlp_shapexx", fixed = TRUE)
  expect_match(two, "nlp_shapexx", fixed = TRUE)
  expect_identical(sub("shapexone", "", one, fixed = TRUE),
                   sub("shapextwo", "", two, fixed = TRUE))

  # And it is DECLARED once, with no per-margin coefficient vectors alongside
  # it. Match the declaration rather than the bare name: `b_shapexx;` also
  # appears where the linear predictor is accumulated
  # (`nlp_shapexx += X_shapexx * b_shapexx;`), so counting occurrences of the
  # name would count a use as a second parameter.
  decl <- grep("vector[K_shapexx] b_shapexx;", lines, fixed = TRUE)
  expect_length(decl, 1L)
  expect_length(grep("vector[K_shapexone] b_shapexone;", lines, fixed = TRUE), 0L)
  expect_length(grep("vector[K_shapextwo] b_shapextwo;", lines, fixed = TRUE), 0L)
})

test_that("the tied model's prior is class 'b' / nlpar, not dpar", {
  p <- as.data.frame(brms::get_prior(tied_formula(), data = sd,
                                     family = binegbin(),
                                     stanvars = binegbin_stanvars()))
  row <- unique(p[p$nlpar == "shapexx", c("class", "dpar", "nlpar")])
  expect_identical(row$class, "b")
  expect_identical(row$dpar,  "")

  # The pre-0.10.0 spelling names nothing in this model, which is exactly why
  # a prior written that way is dropped without complaint.
  expect_equal(sum(p$dpar == "shapex"), 0L)
  expect_equal(sum(p$dpar == "shapexone"), 0L)
})
