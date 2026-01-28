# Tier 4 discrete color scale

Uses the tier4 colorblind-friendly discrete palette (ggthemes colorblind
hex codes).

## Usage

``` r
scale_color_tier4(
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

  Color for NA values.

## Value

A ggplot2 scale.
