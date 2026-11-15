# rbicount(): shape of the result, recycling, the identities y1 = Ns + N1 and
# y2 = Ns + N2, the component moments, the Poisson case, the error for a
# partial set of dispersions, composition with mfd_to_dpars(), and agreement
# with .pair_draw_<family>(). No test here compiles Stan.

NB <- list(
  mu        = 4,
  lambdaone = 1.5,
  lambdatwo = 3,
  shapes    = 2,
  shapexone = 5,
  shapextwo = 0.8
)

test_that("rbicount() returns a data frame of n rows with five named columns", {
  set.seed(1)
  out <- rbicount(7, mu = 2, lambdaone = 1, lambdatwo = 1)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 7L)
  expect_named(out, c("y1", "y2", "Ns", "N1", "N2"))

  out_nb <- do.call(rbicount, c(list(n = 7), NB))
  expect_named(out_nb, c("y1", "y2", "Ns", "N1", "N2"))
  # Integer under both constructors, although rnbinom() returns doubles
  expect_true(all(vapply(out, is.integer, logical(1))))
  expect_true(all(vapply(out_nb, is.integer, logical(1))))

  expect_identical(nrow(rbicount(0, mu = 2, lambdaone = 1, lambdatwo = 1)), 0L)
  # As in rnbinom(), a vector n gives length(n) draws.
  expect_identical(nrow(rbicount(1:4, mu = 2, lambdaone = 1, lambdatwo = 1)), 4L)
})

test_that("rbicount() returns y1 = Ns + N1 and y2 = Ns + N2 on every row", {
  set.seed(2)
  for (out in list(
    rbicount(500, mu = 3, lambdaone = 1, lambdatwo = 2),
    do.call(rbicount, c(list(n = 500), NB))
  )) {
    expect_true(all(out$y1 == out$Ns + out$N1))
    expect_true(all(out$y2 == out$Ns + out$N2))
    expect_true(all(unlist(out) >= 0))
  }
})

test_that("rbicount() recycles each parameter to n", {
  set.seed(3)
  # A rate of zero gives a component of zero, so the recycled pattern can be
  # read off the draws exactly.
  out <- rbicount(
    6,
    mu        = c(0, 50),
    lambdaone = 0,
    lambdatwo = c(0, 0, 50)
  )
  expect_true(all(out$Ns[c(1, 3, 5)] == 0))
  expect_true(all(out$Ns[c(2, 4, 6)] > 0))
  expect_true(all(out$N1 == 0))
  expect_true(all(out$N2[c(1, 2, 4, 5)] == 0))
  expect_true(all(out$N2[c(3, 6)] > 0))

  out_nb <- rbicount(
    6,
    mu        = c(0, 50),
    lambdaone = 0,
    lambdatwo = 0,
    shapes    = c(1, 2, 3),
    shapexone = 1,
    shapextwo = 1
  )
  expect_true(all(out_nb$Ns[c(1, 3, 5)] == 0))
  expect_true(all(out_nb$N1 == 0))
})

# Monte Carlo tolerance: with R draws, the standard error of a sample mean is
# sqrt(var/R), and of a sample variance is about sqrt((mu4 - var^2)/R). The
# tolerance is five standard errors. The fourth central moment of NB2(m, phi)
# is computed from its cumulants.
nb_var <- function(m, phi) m + m^2 / phi
nb_mu4 <- function(m, phi) {
  p  <- phi / (phi + m)
  q  <- 1 - p
  k2 <- phi * q / p^2
  k4 <- phi * q * (1 + 4 * q + q^2) / p^4
  k4 + 3 * k2^2
}

test_that("NB2 components match their means and variances", {
  set.seed(20261115)
  R <- 2e5
  out <- do.call(rbicount, c(list(n = R), NB))
  comps <- list(
    Ns = c(NB$mu,        NB$shapes),
    N1 = c(NB$lambdaone, NB$shapexone),
    N2 = c(NB$lambdatwo, NB$shapextwo)
  )
  for (nm in names(comps)) {
    m   <- comps[[nm]][1]
    phi <- comps[[nm]][2]
    v   <- nb_var(m, phi)
    x   <- out[[nm]]
    expect_lt(abs(mean(x) - m), 5 * sqrt(v / R), label = paste("mean of", nm))
    expect_lt(
      abs(var(x) - v),
      5 * sqrt((nb_mu4(m, phi) - v^2) / R),
      label = paste("variance of", nm)
    )
  }
  # Cov(y1, y2) = Var(Ns)
  v_s <- nb_var(NB$mu, NB$shapes)
  expect_lt(abs(cov(out$y1, out$y2) - v_s), 0.05 * v_s)
})

test_that("Poisson components match their means and variances", {
  set.seed(20261116)
  R <- 2e5
  rates <- c(Ns = 4, N1 = 1.5, N2 = 3)
  out <- rbicount(R, mu = rates[["Ns"]], lambdaone = rates[["N1"]],
                  lambdatwo = rates[["N2"]])
  for (nm in names(rates)) {
    m <- rates[[nm]]
    x <- out[[nm]]
    # Poisson: variance m, fourth central moment m + 3 m^2
    expect_lt(abs(mean(x) - m), 5 * sqrt(m / R), label = paste("mean of", nm))
    expect_lt(abs(var(x) - m), 5 * sqrt((m + 3 * m^2 - m^2) / R),
              label = paste("variance of", nm))
  }
  expect_true(all(out$Ns == round(out$Ns)))
})

test_that("an infinite dispersion gives a Poisson component", {
  set.seed(4)
  R <- 1e5
  out <- rbicount(R, mu = 4, lambdaone = 1, lambdatwo = 1,
                  shapes = Inf, shapexone = Inf, shapextwo = Inf)
  expect_lt(abs(var(out$Ns) - 4), 5 * sqrt((4 + 2 * 16) / R))
})

test_that("rbicount() stops when some dispersions are supplied and not others", {
  expect_error(
    rbicount(5, mu = 1, lambdaone = 1, lambdatwo = 1, shapes = 2),
    "`shapexone`, `shapextwo`"
  )
  expect_error(
    rbicount(5, mu = 1, lambdaone = 1, lambdatwo = 1,
             shapexone = 2, shapextwo = 2),
    "Missing: `shapes`"
  )
  # mfd_to_dpars() with kappas alone writes shapes alone.
  pars <- mfd_to_dpars(M = 5, f = 0.7, kappas = 0.5)
  expect_error(do.call(rbicount, c(list(n = 5), pars)), "all three")
})

test_that("rbicount() composes with mfd_to_dpars()", {
  pars <- mfd_to_dpars(
    M      = 5,
    f      = 0.7,
    delta  = atanh(0.3),
    kappas = 0.5,
    kappax = 0.5
  )
  set.seed(5)
  a <- do.call(rbicount, c(list(n = 50), pars))
  set.seed(5)
  b <- rbicount(
    50,
    mu        = pars$mu,
    lambdaone = pars$lambdaone,
    lambdatwo = pars$lambdatwo,
    shapes    = pars$shapes,
    shapexone = pars$shapexone,
    shapextwo = pars$shapextwo
  )
  expect_identical(a, b)

  # Without overdispersions, mfd_to_dpars() gives the Poisson case.
  set.seed(6)
  p <- do.call(rbicount, c(list(n = 50), mfd_to_dpars(M = 5, f = 0.7)))
  expect_named(p, c("y1", "y2", "Ns", "N1", "N2"))

  # Midpoint of the two expected counts, recovered from the draws
  set.seed(7)
  big <- do.call(rbicount, c(list(n = 1e5), pars))
  expect_lt(abs((mean(big$y1) + mean(big$y2)) / 2 - 5), 0.05)
})

test_that("rbicount() works inside dplyr::mutate() with n() and parameter columns", {
  skip_if_not_installed("dplyr")
  units <- data.frame(treatment = rep(c("a", "b", "c"), each = 30))
  set.seed(8)
  out <- dplyr::mutate(
    units,
    f = c(a = 0.5, b = 0.7, c = 0.9)[treatment],
    as.data.frame(
      mfd_to_dpars(M = 5, f = f, delta = 0, kappas = 0.5, kappax = 0.5)
    ),
    rbicount(dplyr::n(), mu, lambdaone, lambdatwo, shapes, shapexone, shapextwo)
  )
  expect_identical(nrow(out), 90L)
  expect_true(all(c("y1", "y2", "Ns", "N1", "N2", "mu") %in% names(out)))
  expect_true(all(out$y1 == out$Ns + out$N1))

  # Grouped: n() is the size of each group, and the draws are made per group.
  set.seed(9)
  grouped <- units |>
    dplyr::group_by(treatment) |>
    dplyr::mutate(
      rbicount(dplyr::n(), mu = 4, lambdaone = 1, lambdatwo = 1)
    ) |>
    dplyr::ungroup()
  expect_identical(nrow(grouped), 90L)

  # A partially paired design: 20 groups with a random effect on log M, a
  # three-level treatment acting on f and beta, 30 units per cell, and y1
  # recorded on 20% of units.
  set.seed(12)
  u <- rnorm(20, 0, 0.3)
  sim <- expand.grid(
    unit      = 1:30,
    treatment = c("a", "b", "c"),
    group     = 1:20,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      M     = exp(log(5) + u[group]),
      f     = c(a = 0.5, b = 0.7, c = 0.9)[treatment],
      delta = atanh(c(a = 0, b = 0.3, c = -0.3)[treatment]),
      as.data.frame(mfd_to_dpars(M, f, delta, kappas = 0.5, kappax = 0.5)),
      rbicount(dplyr::n(), mu, lambdaone, lambdatwo, shapes, shapexone, shapextwo),
      y1_obs = rbinom(dplyr::n(), 1, 0.2),
      y1     = dplyr::if_else(y1_obs == 1, y1, -1L)
    )
  expect_identical(nrow(sim), 1800L)
  expect_true(all(sim$y1[sim$y1_obs == 0] == -1L))
  expect_true(all(sim$y1[sim$y1_obs == 1] ==
                    (sim$Ns + sim$N1)[sim$y1_obs == 1]))
})

test_that("rbicount() and .pair_draw_<family>() draw the same pairs", {
  # Both functions call .draw_components(), so with one seed the pairs are
  # identical.
  S <- 200
  prep_nb <- make_synthetic_prep(
    dpars = list(
      mu        = rep(NB$mu, S),
      lambdaone = rep(NB$lambdaone, S),
      lambdatwo = rep(NB$lambdatwo, S),
      shapes    = rep(NB$shapes, S),
      shapexone = rep(NB$shapexone, S),
      shapextwo = rep(NB$shapextwo, S)
    ),
    Y     = 0L,
    vint1 = 0L
  )
  set.seed(10)
  pair <- .pair_draw_binegbin(1, prep_nb)
  set.seed(10)
  sim <- do.call(rbicount, c(list(n = S), NB))
  expect_equal(unname(pair[, 1]), sim$y1)
  expect_equal(unname(pair[, 2]), sim$y2)

  prep_p <- make_synthetic_prep(
    dpars = list(
      mu        = rep(NB$mu, S),
      lambdaone = rep(NB$lambdaone, S),
      lambdatwo = rep(NB$lambdatwo, S)
    ),
    Y     = 0L,
    vint1 = 0L
  )
  set.seed(11)
  pair <- .pair_draw_bipois(1, prep_p)
  set.seed(11)
  sim <- rbicount(S, mu = NB$mu, lambdaone = NB$lambdaone,
                  lambdatwo = NB$lambdatwo)
  expect_equal(unname(pair[, 1]), sim$y1)
  expect_equal(unname(pair[, 2]), sim$y2)
})
