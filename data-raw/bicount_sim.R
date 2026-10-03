# Generates `bicount_sim` and `bicount_sim_latent`, a simulated partially paired
# dataset for the examples and the methods paper.
#
#   Rscript data-raw/bicount_sim.R      # from the package root
#
# The design is generic and is not tuned to any real dataset:
#   - three sites, 180 units each;
#   - source 2 recorded on every unit;
#   - source 1 recorded on exactly 36 units per site (20%), drawn at random within
#     site, so the pairing is missing completely at random by construction;
#   - one pair of counts per unit from the joint model with negative-binomial
#     components, with site effects on congruence f and on bias beta = tanh(delta).
#
# The generating values are stated once, below, and repeated in R/data.R.

SEED <- 20261002L
set.seed(SEED)

n_per_site <- 180L
pair_frac  <- 0.2

truth <- data.frame(
  site    = factor(c("A", "B", "C")),
  M       = 5,
  f       = c(0.5, 0.7, 0.9),
  beta    = c(0, 0.3, -0.3),
  kappas  = 0.5,
  kappax1 = 0.5,
  kappax2 = 0.5
)
truth$delta   <- atanh(truth$beta)
truth$mu      <- truth$M * truth$f
truth$lambda1 <- truth$M * (1 - truth$f) * (1 + truth$beta)
truth$lambda2 <- truth$M * (1 - truth$f) * (1 - truth$beta)

## Generated site by site, in a fixed order, from the single seed above.
sites <- lapply(seq_len(nrow(truth)), function(j) {
  tr <- truth[j, ]
  Ns <- stats::rnbinom(n_per_site, size = 1 / tr$kappas^2,  mu = tr$mu)
  N1 <- stats::rnbinom(n_per_site, size = 1 / tr$kappax1^2, mu = tr$lambda1)
  N2 <- stats::rnbinom(n_per_site, size = 1 / tr$kappax2^2, mu = tr$lambda2)
  obs <- integer(n_per_site)
  obs[sample.int(n_per_site, round(pair_frac * n_per_site))] <- 1L
  data.frame(
    site   = tr$site,
    unit   = seq_len(n_per_site),
    y1     = ifelse(obs == 1L, as.integer(Ns + N1), NA_integer_),
    y2     = as.integer(Ns + N2),
    y1_obs = obs,
    Ns     = as.integer(Ns),
    N1     = as.integer(N1),
    N2     = as.integer(N2)
  )
})
all <- do.call(rbind, sites)

bicount_sim        <- all[, c("site", "unit", "y1", "y2", "y1_obs")]
bicount_sim_latent <- all[, c("site", "unit", "Ns", "N1", "N2")]

## ---- checks -------------------------------------------------------------------
lat_y1 <- bicount_sim_latent$Ns + bicount_sim_latent$N1
lat_y2 <- bicount_sim_latent$Ns + bicount_sim_latent$N2
frac   <- tapply(bicount_sim$y1_obs, bicount_sim$site, mean)
stopifnot(
  "wrong number of rows" = nrow(bicount_sim) == 3L * n_per_site,
  "the two objects are not row-aligned" =
    identical(bicount_sim[, c("site", "unit")], bicount_sim_latent[, c("site", "unit")]),
  "y1 is not Ns + N1 on every recorded unit" =
    all(bicount_sim$y1[bicount_sim$y1_obs == 1L] == lat_y1[bicount_sim$y1_obs == 1L]),
  "y2 is not Ns + N2 on every unit" = all(bicount_sim$y2 == lat_y2),
  "y1 is not NA exactly where y1_obs == 0" =
    all(is.na(bicount_sim$y1) == (bicount_sim$y1_obs == 0L)),
  "pairing fraction more than 2 percentage points from its target in a site" =
    all(abs(frac - pair_frac) <= 0.02)
)

usethis::use_data(bicount_sim, bicount_sim_latent, overwrite = TRUE,
                  compress = "xz", version = 2)
