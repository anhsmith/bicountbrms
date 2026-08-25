# Convert (M, f, delta) coordinates to native binegbin/bipois dpars

Maps the interpretable coordinates – overall level \`M\`, congruence
\`f\`, and source bias \`delta\` – onto the rate dpars every family in
this package takes (\`mu\`, \`lambdaone\`, \`lambdatwo\`), optionally
converting SD-scale dispersions to the
\`shapes\`/\`shapexone\`/\`shapextwo\` dpars.

The three rates are common to \[bipois()\], \[bipois_partialobs()\],
\[binegbin()\] and \[binegbin_partialobs()\], so this direction serves
all four constructors. The dispersions are where they differ: the
Poisson families have none, and both negative-binomial constructors have
one per margin.

Everything this returns is named for a dpar a shipping family accepts,
so the output can go straight into a \`brm()\` call or a simulation.
That is why \`kappax\` writes \`shapexone\` and \`shapextwo\` at a
common value rather than a single \`shapex\`: no family has had a
\`shapex\` dpar since 0.10.0.

\[binegbin_dpars_to_mfd()\] is the inverse. It still accepts \`shapex\`,
because its input is a stored fit and pre-0.10.0 \`binegbin\` fits carry
that name.

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

  Optional SD-scale dispersions. \`0\` is the Poisson limit. \`kappas\`
  is the shared component's, and the returned list gains \`shapes\` (\`=
  1/kappa^2\`, so \`kappa = 0\` gives \`Inf\`). \`kappax\` is the
  shorthand for a single excess dispersion governing \*both\* margins:
  supply it and the returned list gains \`shapexone\` and \`shapextwo\`
  at that common value, which is the symmetric model \[binegbin()\]
  reaches by tying the two with \`nlf()\`. Omit both for \[bipois()\]
  and \[bipois_partialobs()\], which have no dispersion parameters.

- kappaxone, kappaxtwo:

  Optional per-margin SD-scale excess dispersions, for the general case
  in which the two margins are free to differ. If supplied, the returned
  list gains \`shapexone\`/\`shapextwo\`. Mutually exclusive with
  \`kappax\`, which writes the same two slots.

## Value

A named list of \`mu\`, \`lambdaone\`, \`lambdatwo\`, plus \`shapes\`
when \`kappas\` is supplied, and \`shapexone\`/\`shapextwo\` when either
\`kappax\` or \`kappaxone\`/\`kappaxtwo\` are. Every name is a dpar of a
shipping family.

## Details

Arguments are recycled to a common length, so this vectorises over
posterior draws.

\*\*Boundary behaviour.\*\* At \`f = 1\` both excess rates are exactly
\`0\` regardless of \`delta\` – the bias becomes unidentifiable, which
\[binegbin_dpars_to_mfd()\] reports back as \`NA\`. This direction is
always well defined; only the inverse degenerates.

\*\*The two ways to name the excess dispersion.\*\* Both
negative-binomial constructors carry the pair
\`shapexone\`/\`shapextwo\`, so that is what this function writes.
Supply \`kappaxone\`/\`kappaxtwo\` to give the two margins different
values, or \`kappax\` to give them the same one – the symmetric model,
term for term the pre-0.8.0 likelihood, which a fit reaches by routing
both dpars through one non-linear parameter. The two spellings are
mutually exclusive in a single call because they write the same two
slots.

Before 0.10.0, \`kappax\` returned a dpar named \`shapex\` and
\`kappaxone\`/\`kappaxtwo\` returned the pair, because two different
families wanted two different things. There is now one family and one
dpar set, so both spellings produce it.

\*\*Using these coordinates with the Poisson families.\*\* \[bipois()\]
and \[bipois_partialobs()\] take the same three rates and no dispersion,
so call this with \`M\`, \`f\` and \`delta\` alone and pass the result
straight through. There is no \`kappa\` to supply: the Poisson case is
not a dispersion set to a particular value but the absence of the
parameter, which is precisely why fitting it wants its own family rather
than \[binegbin()\] with \`kappa\` driven to \`0\`.

## See also

\[binegbin_dpars_to_mfd()\], \[bipois()\], \[bipois_partialobs()\],
\[binegbin()\], \[binegbin_partialobs()\]

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

# One excess dispersion for both margins, and two different ones
binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4, kappax = 0.9)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 3.96
#> 
#> $lambdatwo
#> [1] 3.96
#> 
#> $shapes
#> [1] 6.25
#> 
#> $shapexone
#> [1] 1.234568
#> 
#> $shapextwo
#> [1] 1.234568
#> 
binegbin_mfd_to_dpars(M = 12, f = 0.67, kappas = 0.4,
                      kappaxone = 0.9, kappaxtwo = 0.3)
#> $mu
#> [1] 8.04
#> 
#> $lambdaone
#> [1] 3.96
#> 
#> $lambdatwo
#> [1] 3.96
#> 
#> $shapes
#> [1] 6.25
#> 
#> $shapexone
#> [1] 1.234568
#> 
#> $shapextwo
#> [1] 11.11111
#> 

# To FIT in these coordinates, pass them through a non-linear formula
# (every dpar is log-linked, so the link supplies the exp()):
#   bf(y1 | vint(y2) ~ 1, nl = TRUE) +
#     nlf(mu        ~ eta + log_inv_logit(con)) +
#     nlf(lambdaone ~ log(2) + eta + log_inv_logit(-con) +
#                     log_inv_logit(2 * methd)) +
#     nlf(lambdatwo ~ log(2) + eta + log_inv_logit(-con) +
#                     log_inv_logit(-2 * methd)) +
#     lf(eta ~ 1, con ~ 1, methd ~ 1)
# where eta = log M, con = logit f, methd = delta.
```
