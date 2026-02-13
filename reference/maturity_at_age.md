# Maturity-at-age

Build a maturity-at-age vector from a maturity specification.

## Usage

``` r
maturity_at_age(ages, maturity)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- maturity:

  A named list specifying maturity. Required element `type`. For
  `type="vector"`, provide `m`. For `type="logistic"`, provide `a50` and
  `delta`.

## Value

Numeric vector of length `length(ages)` with values in \\\[0,1\]\\.

## Details

This is a user-facing wrapper around an internal helper used by
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
It is useful for sensitivity analyses where the user wants to modify the
maturity specification of an existing
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
object.

It returns maturity bounded to \\\[0,1\]\\.
