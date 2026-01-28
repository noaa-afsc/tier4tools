# Construct SPR input object

Validates and expand inputs needed for per-recruit SPR/YPR calculations
for one or multiple species.

## Usage

``` r
spr_input(
  ages,
  species,
  rec_prop = NULL,
  R0 = 1,
  use_plus_group = TRUE,
  scale_selex = TRUE
)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- species:

  A named or unnamed list for single species specification, or a named
  list of species for multispecies specifications. See details.

- rec_prop:

  Optional numeric vector of recruitment proportions for multispecies
  analyses. Must be the same length as the number of species and sum
  to 1. Defaults to equal proportions.

- R0:

  Numeric scalar recruitment scaling used for per-recruit calculations
  (default=1, or equal proportions if inputting multiple species). If
  defining R0 for multiple species, R0 proportions must be a vector
  equal to the number of species and sum to 1.

- use_plus_group:

  Logical. If `TRUE` (default), the terminal age will be treated as a
  plus group in pop dy equations.

- scale_selex:

  Logical. If `TRUE` (default), scales selectivity so max equals 1 when
  possible.

## Value

A named list with expanded, validated inputs for running SPR functions:

- `ages`, `R0`, `use_plus_group`

- `species`: named list of per-species expanded vectors and original
  specifications

- `rec_prop`: recruitment proportions used for multispecies calculations

## Details

For a single-species analysis, provide a single species specification
via `species` (can be named or unnamed). For multispecies analyses,
provide a named list of species specifications, where each element
defines life history and fishery inputs for one species.

Each species specification must include:

- `len_at_age`: list with `type="vb"` (Linf, k, t0) or `type="vector"`
  (L).

- `wt_at_age`: list with `type="wl"` (alpha, beta) or `type="vector"`
  (W).

- `maturity`: list with `type="logistic"` (a50, delta) or
  `type="vector"` (m).

- `selectivity`: list with `type="logistic"` (a50, delta),
  `type="knife_edge"` (a_c), or `type="vector"` (s).

- `M`: numeric scalar or vector of length `length(ages)`.

Rules for size-at-age:

- If `wt_at_age$type = "vector"`, the provided weight-at-age vector is
  used directly in the SPR model. In this case, `len_at_age` is not used
  to determine size-at-age for the SPR calculations, and a message is
  emitted to make this explicit.

- If `len_at_age$type = "vb"`, then `wt_at_age$type` must be `"wl"` to
  ensure consistency between growth and weight conversion.

## Examples
