# Fishery selectivity-at-age

Build a selectivity-at-age vector from a selectivity specification.

## Usage

``` r
selectivity_at_age(ages, selectivity, scale_max = TRUE)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- selectivity:

  A named list specifying selectivity. Required element `type`. For
  `type="vector"`, provide `s`. For `type="knife_edge"`, provide `a_c`.
  For `type="logistic"`, provide `a50` and `delta`.

- scale_max:

  Logical. If `TRUE`, scales selectivity so that `max(s_a) = 1` when
  possible.

## Value

Numeric vector of length `length(ages)` with values in \\\[0,1\]\\. When
`selectivity$type == "fleets"`, the returned vector is the effective
selectivity-at-age.

## Details

This is a user-facing wrapper around an internal helper used by
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
It is useful for sensitivity analyses where the user wants to modify the
fishery selectivity specification of an existing
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
object.

It returns selectivity bounded to \\\[0,1\]\\ and can optionally scale
so that `max(s_a) = 1`.
