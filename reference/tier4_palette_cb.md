# Tier 4 default colorblind-friendly palette (ggthemes "colorblind" hex codes)

Discrete palette intended for species, scenarios, or other categorical
aesthetics.

## Usage

``` r
tier4_palette_cb(n, reverse = FALSE, recycle = TRUE, warn = TRUE)
```

## Arguments

- n:

  Number of colors requested.

- reverse:

  Logical. If TRUE, reverse the palette order.

- recycle:

  Logical. If TRUE and n exceeds palette length, recycle colors (default
  TRUE).

- warn:

  Logical. If TRUE, warn when recycling is needed (default TRUE).

## Value

Character vector of hex color codes.
