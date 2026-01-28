# Run per-recruit SPR and YPR calculations

Computes spawning potential ratio (SPR) and yield per recruit (YPR)
across a range of fishing mortality values for a single stock or
multispecies complex defined by
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
The function returns species-specific and combined (total) SPR(F) and
YPR(F) curves, SPR-based fishing mortality reference points (e.g., F40,
F35), and optional decomposition diagnostics.

## Usage

``` r
run_spr(
  x,
  F_vec = seq(0, 0.4, by = 0.001),
  spr_targets = c(0.4, 0.35),
  multispecies_constraint = c("none", "all_species"),
  diagnostics = TRUE
)
```

## Arguments

- x:

  An object returned by
  [`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
  defining ages and species-specific life history schedules.

- F_vec:

  Numeric vector of fishing mortality values (F) to evaluate. Must be
  finite and \>= 0.

- spr_targets:

  Numeric vector of SPR targets in (0,1), for example `c(0.40, 0.35)`.

- multispecies_constraint:

  Character string specifying whether to apply a multispecies constraint
  when selecting combined reference points. Options are:

  `"none"`

  :   Unconstrained reference points are based on the combined SPR curve
      only.

  `"all_species"`

  :   Constrained combined reference points, all species must meet the
      target SPR.

- diagnostics:

  Logical. If TRUE, returns plus-group diagnostics and age-specific
  decompositions.

## Value

A named list with the following elements:

- `F_vec`:

  Fishing mortality grid used for calculations.

- `spr_targets`:

  SPR targets requested.

- `multispecies_constraint`:

  Constraint setting used.

- `SPR_total`:

  Combined SPR(F) across species (length = `length(F_vec)`).

- `YPR_total`:

  Combined YPR(F) across species (length = `length(F_vec)`).

- `SPR_by_species`:

  Named list of species-specific SPR(F) vectors.

- `YPR_by_species`:

  Named list of species-specific YPR(F) vectors.

- `F_spr_total`:

  Named numeric vector of combined reference points (nearest-grid
  match), for example `F40_total` and `F35_total`.

- `F_spr_total_constrained`:

  Named numeric vector of constrained combined reference points when
  `multispecies_constraint="all_species"`; otherwise `NULL`.

- `F_spr_by_species`:

  Named list of species-specific reference points for each SPR target.

- `diagnostics`:

  A list with elements:

  - `Ftarget`: the fishing mortality value used for diagnostics
    (defaults to the first target, and uses constrained `Ftarget` when
    available),

  - `plus_share`: plus-group share table (or `NULL` if
    `diagnostics=FALSE`),

  - `decomp`: age-specific decomposition data frame (or `NULL` if
    `diagnostics=FALSE`).

## Details

**SPR and YPR.** For each species, unfished spawning biomass per recruit
(SBPR0) and fished SBPR(F) are computed from survivorship-at-age,
weight-at-age, and maturity-at-age. SPR is defined as
`SPR(F) = SBPR(F) / SBPR(0)`. YPR is computed using the Baranov catch
equation in biomass units and the age-specific fishery selectivity.

**Multispecies complexes.** For multispecies inputs, each species is
evaluated separately using its own life history schedules and
recruitment scaling (i.e., recruitment proportions in
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)).
Combined SBPR and yield are computed by summing across species.

**SPR reference points.** Fishing mortality reference points are
reported for:

- the combined complex (`F_spr_total`), based on matching the combined
  SPR curve to each value in `spr_targets` using a nearest-grid
  approach, and

- each species individually (`F_spr_by_species`).

**Multispecies constraint.** When
`multispecies_constraint = "all_species"`, a constrained reference point
(`F_spr_total_constrained`) is calculated for each SPR target as the max
F in `F_vec` such that *both* the combined SPR and each species-specific
SPR meet or exceed the target. If no value satisfies the constraint,
`NA` is returned for that target.

**Diagnostics and decomposition outputs.** When `diagnostics = TRUE`,
the output includes:

- `diagnostics$plus_share`: the fraction of total SBPR contributed by
  the terminal age (often a plus group) at F = 0 and at the diagnostic
  reference F `Ftarget`, and

- `diagnostics$decomp`: an age-specific df containing survivorship,
  maturity, selectivity, weight-at-age, M-at-age, SBPR contributions,
  and SBPR removed (F0 minus Ftarget) for each species at F = 0 and F =
  Ftarget.

These diagnostics can be visualized with
[`plot_spr_decomp()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_decomp.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume inp is a spr_input() object (single species or multispecies).
# Run per-recruit calculations and request decomposition diagnostics.
spr_out <- run_spr(inp, diagnostics = TRUE,
                  multispecies_constraint = "all_species")

# reference points
spr_out$F_spr_total
spr_out$F_spr_total_constrained
spr_out$F_spr_by_species

# Plot SPR(F) curves, complex-level SPR and species-specific SPRs
# plot_spr_curves(spr_out, which = "total")
# plot_spr_curves(spr_out, which = "species")

# Decomposition diagnostics
plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
plots$contrib
plots$removed
plots$survivorship

# Plus-group diagnostic table (fraction of SBPR in terminal age)
spr_out$diagnostics$plus_share
} # }
```
