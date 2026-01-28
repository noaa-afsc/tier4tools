# Length-at-age

Convert a user-specified length-at-age specification into an
age-specific vector spanning the full age range `ages = a_min:a_max`.

## Usage

``` r
calc_len_at_age(ages, len_at_age)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- len_at_age:

  A named list specifying length-at-age. Required element `type`. For
  `type="vb"`, provide `Linf`, `k`, and `t0`. For `type="vector"`,
  provide `L`.

## Value

Numeric vector `L_a` with `length(ages)`.

## Details

Supported inputs:

- `type = "vb"`: von Bertalanffy growth with parameters `Linf`, `k`,
  `t0`.

- `type = "vector"`: external length-at-age vector `L` with
  `length(L) == length(ages)`.

Von Bertalanffy form: \$\$L_a = L\_\infty \left(1 - \exp\left\[-k(a -
t_0)\right\]\right)\$\$
