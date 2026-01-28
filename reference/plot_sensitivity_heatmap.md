# Plot sensitivity heatmap for SPR outputs

Heatmap based on
[`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md)
output for a two-parameter sweep with the base values overlayed. Uses
viridis by default when available.

## Usage

``` r
plot_sensitivity_heatmap(
  sens_df,
  x,
  y,
  fill,
  facet = NULL,
  title = NULL,
  fill_scale = c("viridis", "base"),
  base_marker = TRUE,
  base_point_size = 5,
  base_point_color = "#D55E00",
  base_point_stroke = 1.4,
  base_label = FALSE
)
```

## Arguments

- sens_df:

  Data frame returned by
  [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md).

- x:

  Name of the column to use for the x-axis (character).

- y:

  Name of the column to use for the y-axis (character).

- fill:

  Name of the column to use for tile fill (character).

- facet:

  Optional column name to facet by (character).

- title:

  Optional plot title.

- fill_scale:

  Character, either `"base"` or `"viridis"` (default).

- base_marker:

  Logical. If TRUE, overlays a base-case marker when base values are
  available.

- base_point_size:

  Numeric, size of base case mark (X).

- base_point_color:

  Outline color for base case mark, defaults to a colorblind-friendly
  orange.

- base_point_stroke:

  Stroke width for base case mark.

- base_label:

  Logical. If TRUE, labels the base marker with "Base".

## Value

A ggplot object.

## Details

This function is intended for visualization of structured parameter
sweeps where `x` and `y` are numeric columns in the sensitivity table
and `fill` is a numeric output (for example `"F40_total"`,
`"F40_total_constrained"`, `"YPR_at_F40"`, or
`"plus_share_RE_Ftarget"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Heatmap of unconstrained F40 across (M, a50)
plot_sensitivity_heatmap(out, x = "M", y = "maturity_a50", fill = "F40_total")
} # }
```
