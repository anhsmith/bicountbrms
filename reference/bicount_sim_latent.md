# Latent components of the simulated partially paired counts

The three latent components behind each pair of counts in
[bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md):
the count recorded by both sources and the counts recorded by each
source alone. Kept apart from
[bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md)
because they would not be observed in practice; they give the true value
of every count, including the source-1 counts that were not recorded,
against which estimates and predictions can be checked.

## Usage

``` r
bicount_sim_latent
```

## Format

A data frame with 540 rows, aligned row for row with
[bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md),
and 5 columns:

- site, unit:

  As in
  [bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md).

- Ns:

  Integer, the count recorded by both sources.

- N1:

  Integer, the count recorded by source 1 only.

- N2:

  Integer, the count recorded by source 2 only.

## Details

`y1 = Ns + N1` on every unit, recorded or not, and `y2 = Ns + N2`.
Generating values and seed are given in
[bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md).

## See also

[bicount_sim](https://anhsmith.github.io/bicountbrms/reference/bicount_sim.md)

## Examples

``` r
# The true source-1 total over the units on which source 1 was not recorded
un <- bicount_sim$y1_obs == 0
tapply(bicount_sim_latent$Ns[un] + bicount_sim_latent$N1[un],
       bicount_sim$site[un], sum)
#>   A   B   C 
#> 661 805 718 
```
