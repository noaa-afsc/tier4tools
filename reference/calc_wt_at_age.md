# Weight-at-age

Convert a user-specified weight-at-age specification into an
age-specific vector spanning the full age range `ages = a_min:a_max`.

## Usage

``` r
calc_wt_at_age(ages, wt_at_age, len_at_age = NULL)
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

Numeric vector `w_a` with `length(ages)`.

## Details

Supported inputs:

- `type = "wl"`: derive weight-at-age from length-at-age using `alpha`
  and `beta`.

- `type = "vector"`: external weight-at-age vector `W` with
  `length(W) == length(ages)`.

Weight-length form: \$\$w_a = \alpha L_a^\beta\$\$
