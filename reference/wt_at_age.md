# Weight-at-age

Build a weight-at-age vector from a weight specification.

## Usage

``` r
wt_at_age(ages, wt_at_age, len_at_age = NULL)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- wt_at_age:

  A named list specifying weight-at-age. Required element `type`. For
  `type="wl"`, provide `alpha` and `beta`. For `type="vector"`, provide
  `W`.

- len_at_age:

  Optional numeric vector `L_a` with `length(ages)`. Required when
  `wt_at_age$type = "wl"`.

## Value

Numeric vector of length `length(ages)`.

## Details

This is a user-facing wrapper around an internal helper used by
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
It is useful for sensitivity analyses where the user wants to modify the
weight specification of an existing
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
object.

Use `len_at_age` when `wt_at_age$type == "wl"`.
