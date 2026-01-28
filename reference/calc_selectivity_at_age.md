# Fishery selectivity-at-age

Convert a user-specified selectivity specification into an age-specific
vector spanning the full age range `ages = a_min:a_max`. Selectivity is
bounded to \\\[0, 1\]\\ and may optionally be scaled so that
`max(s_a) = 1`.

## Usage

``` r
calc_selectivity_at_age(ages, selectivity, scale_max = TRUE)
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

Numeric vector `selectivity_at_age` with `length(ages)`.

## Details

Supported inputs:

- `type = "vector"`: external selectivity vector `s` with
  `length(s) == length(ages)`.

- `type = "knife_edge"`: step function with age-at-full-selection `a_c`.

- `type = "logistic"`: logistic curve with parameters `a50` and `delta`.

Knife-edge form: \$\$s_a = 0 \text{ for } a \< a_c,\quad s_a = 1 \text{
for } a \ge a_c\$\$

Logistic form: \$\$s_a = \frac{1}{1 + \exp\left(\frac{a -
a\_{50}}{\delta}\right)}\$\$
