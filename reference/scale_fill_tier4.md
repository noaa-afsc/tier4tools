# Tier 4 discrete fill scale

Uses the tier4 colorblind-friendly discrete palette (ggthemes colorblind
hex codes).

## Usage

``` r
scale_fill_tier4(
  ...,
  reverse = FALSE,
  recycle = TRUE,
  warn = TRUE,
  na.value = "grey50"
)
```

## Arguments

- ..., reverse, recycle, warn:

  Passed to internal palette helper.

- na.value:

  Fill color for NA values.

## Value

A ggplot2 scale.
