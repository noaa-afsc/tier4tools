# Maturity-at-age

Convert a user-specified maturity specification into an age-specific
vector spanning the full age range `ages = a_min:a_max`.

## Usage

``` r
calc_maturity_at_age(ages, maturity)
```

## Arguments

- ages:

  Numeric vector of ages (e.g., `a_min:a_max`).

- maturity:

  A named list specifying maturity. Required element `type`. For
  `type="vector"`, provide `m`. For `type="logistic"`, provide `a50` and
  `delta`.

## Value

Numeric vector `maturity_at_age` with `length(ages)`.

## Details

Supported inputs:

- `type = "vector"`: external maturity vector `m` with
  `length(m) == length(ages)`.

- `type = "logistic"`: logistic curve with parameters `a50` and `delta`.

Logistic form: \$\$m_a = \frac{1}{1 + \exp\left(\frac{a -
a\_{50}}{\delta}\right)}\$\$
