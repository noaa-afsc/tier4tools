# Length-at-age

Build a length-at-age vector from a length specification.

## Usage

``` r
len_at_age(ages, len_at_age)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- len_at_age:

  A named list specifying length-at-age. Required element `type`. For
  `type="vb"`, provide `Linf`, `k`, and `t0`. For `type="vector"`,
  provide `L`.

## Value

Numeric vector of length `length(ages)`.

## Details

This is a user-facing wrapper around an internal helper used by
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
It is useful for sensitivity analyses where the user wants to modify the
length-at-age specification of an existing
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
object.
