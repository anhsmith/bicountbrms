# Convert (M, f, delta) coordinates to native binegbin/bipois dpars

Maps the interpretable coordinates – overall level \`M\`, congruence
\`f\`, and source bias \`delta\` – onto the rate dpars that
\[binegbin()\] and \[bipois()\] actually take (\`mu\`, \`lambdaone\`,
\`lambdatwo\`), optionally converting SD-scale dispersions
\`kappas\`/\`kappax\` to the \`shapes\`/\`shapex\` dpars.

\[binegbin_dpars_to_mfd()\] is the exact inverse.

## Usage

``` r
binegbin_mfd_to_dpars(
  M,
  f,
  delta = 0,
  kappas = NULL,
  kappax = NULL,
  kappaxone = NULL,
  kappaxtwo = NULL
)
```

## Arguments

- M:

  Overall level: \`mu + (lambdaone + lambdatwo)/2\`. Non-negative.

- f:

  Congruence, the share of \`M\` that both sources saw: \`mu / M\`. In
  \`\[0, 1\]\`. \`f = 1\` means perfect agreement (both excesses
  vanish); \`f = 0\` means no shared component at all.

- delta:

  Source bias on the log-ratio scale, \`0.5 \*
  log(lambdaone/lambdatwo)\`. \`0\` is unbiased. \`+/-Inf\` is permitted
  and gives the limit where one excess rate is zero.

- kappas, kappax:

  Optional SD-scale dispersions for the shared and excess components.
  \`0\` is the Poisson limit. If supplied, the returned list gains
  \`shapes\`/\`shapex\` (\`= 1/kappa^2\`, so \`kappa = 0\` gives
  \`Inf\`). \`kappax\` is the single excess dispersion of
  \[binegbin()\]/\[bipois()\].

- kappaxone, kappaxtwo:

  Optional per-margin SD-scale excess dispersions, for
  \[binegbin_joint()\], which carries one per margin rather than one
  shared. If supplied, the returned list gains
  \`shapexone\`/\`shapextwo\`. Mutually exclusive with \`kappax\`.

## Value

A named list of \`mu\`, \`lambdaone\`, \`lambdatwo\`, plus \`shapes\`
when \`kappas\` is supplied, \`shapex\` when \`kappax\` is, and
\`shapexone\`/\`shapextwo\` when \`kappaxone\`/\`kappaxtwo\` are.

## Details

Arguments are recycled to a common length, so this vectorises over
posterior draws.

\*\*Boundary behaviour.\*\* At \`f = 1\` both excess rates are exactly
\`0\` regardless of \`delta\` – the bias becomes unidentifiable, which
\[binegbin_dpars_to_mfd()\] reports back as \`NA\`. This direction is
always well defined; only the inverse degenerates.

\*\*Using these coordinates with \[binegbin_joint()\].\*\* The three
rates carry over unchanged. Since 0.8.0 that family's excess dispersion
is the pair \`shapexone\`/\`shapextwo\` rather than a single \`shapex\`,
so supply \`kappaxone\`/\`kappaxtwo\` instead of \`kappax\` and the
returned list is named to match its dpars. \`kappax\` and the pair are
mutually exclusive in one call – they are two spellings of the same
quantity for different families, and returning both would leave the
caller to guess which their family wants.

For the symmetric model (one excess dispersion, term for term the
pre-0.8.0 likelihood) pass the same value twice: \`kappaxone = k,
kappaxtwo = k\`.

## See also

\[binegbin_dpars_to_mfd()\], \[binegbin()\], \[bipois()\]

## Examples

``` r
# A moderately congruent pair, source 1 running high
binegbin_mfd_to_dpars(M = 12, f = 0.67, delta = 0.2)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 4.741606
#> 
#> $lambdatwo
#> [1] 3.178394
#> 

# Perfect congruence: both excesses vanish
binegbin_mfd_to_dpars(M = 12, f = 1, delta = 0.5)
#> $mu
#> [1] 12
#> 
#> $lambdaone
#> [1] 0
#> 
#> $lambdatwo
#> [1] 0
#> 

# To FIT in these coordinates, pass them through a non-linear formula
# (all five dpars are log-linked, so the link supplies the exp()):
#   bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#     nlf(mu        ~ eta + log_inv_logit(con)) +
#     nlf(lambdaone ~ log(2) + eta + log_inv_logit(-con) +
#                     log_inv_logit(2 * methd)) +
#     nlf(lambdatwo ~ log(2) + eta + log_inv_logit(-con) +
#                     log_inv_logit(-2 * methd)) +
#     lf(eta ~ 1, con ~ 1, methd ~ 1)
# where eta = log M, con = logit f, methd = delta.
```
