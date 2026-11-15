# tests/testthat/test-stancode-shape.R
#
# The generated Stan call of each constructor. Each family has one Stan lpmf
# that takes the observation indicator, and the two constructors differ only
# in `vars`:
#
#   binegbin()             vars = c("vint1[n]", "1")
#   binegbin_partialobs()  vars = c("vint1[n]", "vint2[n]")
#
# brms pastes the literal "1" into the generated call, so a fully paired model
# calls the same lpmf with y1_obs fixed at 1, and no Stan function overloading
# is required. This behaviour of brms:::stan_log_lik_custom() is not
# documented, so these tests check the generated call directly. stancode() and
# standata() do not compile, so the file requires brms but no Stan toolchain.

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
# The plain constructors pass the observation indicator as a literal
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
  # A model from binegbin() has no second integer column in its Stan data.
  sdta <- brms::standata(
    brms::bf(y1 | vint(y2) ~ 1, lambdaone ~ 1, lambdatwo ~ 1,
             shapes ~ 1, shapexone ~ 1, shapextwo ~ 1),
    data = sd, family = binegbin(), stanvars = binegbin_stanvars()
  )
  expect_true("vint1" %in% names(sdta))
  expect_false("vint2" %in% names(sdta))
})

# ---------------------------------------------------------------------------
# The _partialobs constructors pass the observation indicator as data
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
  # brms builds the Stan function name and the parameter block from the name,
  # dpars, links, and bounds, so these must match; only `vars` differs.
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
  # The _partialobs_stanvars() functions return the same Stan code as the
  # plain versions.
  expect_identical(binegbin_stanvars(), binegbin_partialobs_stanvars())
  expect_identical(bipois_stanvars(),   bipois_partialobs_stanvars())
})

test_that("both negative-binomial constructors declare the full six dpars", {
  want <- c("mu", "lambdaone", "lambdatwo", "shapes", "shapexone", "shapextwo")
  expect_identical(binegbin()$dpars, want)
  expect_identical(binegbin_partialobs()$dpars, want)
  # No constructor declares `shapex`.
  expect_false("shapex" %in% binegbin()$dpars)
})

# ---------------------------------------------------------------------------
# One dispersion for both exclusive components
# ---------------------------------------------------------------------------
#
# ?binegbin documents a formula that ties shapexone and shapextwo through one
# non-linear parameter:
#
#   nlf(shapexone ~ shapexx), nlf(shapextwo ~ shapexx), shapexx ~ 1
#
# The tests check that the generated Stan code assigns both dpars from one
# parameter, and that the prior on that parameter has class = "b" and
# nlpar = "shapexx". A prior with class = "Intercept" and dpar = "shapex"
# names no parameter, and brms stops with an error.

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

  # Both dpars are assigned from the same non-linear parameter.
  expect_match(one, "nlp_shapexx", fixed = TRUE)
  expect_match(two, "nlp_shapexx", fixed = TRUE)
  expect_identical(sub("shapexone", "", one, fixed = TRUE),
                   sub("shapextwo", "", two, fixed = TRUE))

  # b_shapexx is declared once, and no coefficient vectors are declared for
  # shapexone or shapextwo. The pattern matches the declaration, because
  # `b_shapexx;` also appears in `nlp_shapexx += X_shapexx * b_shapexx;`.
  decl <- grep("vector[K_shapexx] b_shapexx;", lines, fixed = TRUE)
  expect_length(decl, 1L)
  expect_length(grep("vector[K_shapexone] b_shapexone;", lines, fixed = TRUE), 0L)
  expect_length(grep("vector[K_shapextwo] b_shapextwo;", lines, fixed = TRUE), 0L)
})

test_that("the prior in the tied model is class 'b' / nlpar, not dpar", {
  p <- as.data.frame(brms::get_prior(tied_formula(), data = sd,
                                     family = binegbin(),
                                     stanvars = binegbin_stanvars()))
  row <- unique(p[p$nlpar == "shapexx", c("class", "dpar", "nlpar")])
  expect_identical(row$class, "b")
  expect_identical(row$dpar,  "")

  # A prior on dpar = "shapex" names no parameter of this model, and brms
  # stops with an error rather than fitting with a flat prior.
  expect_equal(sum(p$dpar == "shapex"), 0L)
  expect_equal(sum(p$dpar == "shapexone"), 0L)
  expect_error(
    brms::stancode(tied_formula(), data = sd, family = binegbin(),
                   stanvars = binegbin_stanvars(),
                   prior = brms::prior(normal(0, 1.5), class = "Intercept",
                                       dpar = "shapex")),
    "do not correspond to any model parameter"
  )
})
