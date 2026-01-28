# Arrange SPR decomposition plots into report-ready grid

Combined list output from
[`plot_spr_decomp()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_decomp.md)
into a single multipanel figure using patchwork. Because patchwork is in
Suggests, this function will only run when it is installed.

## Usage

``` r
grid_spr_decomp(
  plots,
  order = NULL,
  ncol = 1,
  nrow = NULL,
  guides = c("collect", "keep"),
  title = NULL,
  caption = NULL
)
```

## Arguments

- plots:

  A named list of ggplot objects returned by
  [`plot_spr_decomp()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_decomp.md).

- order:

  Character vector giving the plot names (in `plots`) in the order to
  arrange. Defaults to
  `c("contrib","removed","removed_prop","survivorship")` when present,
  otherwise `names(plots)`.

- ncol:

  Number of columns in the grid.

- nrow:

  Optional number of rows. If `NULL`, patchwork chooses automatically.

- guides:

  Passed to `patchwork::plot_layout(guides = ...)`. Default `"collect"`
  collects legends when possible.

- title:

  Optional figure title (character). Added with
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html).

- caption:

  Optional figure caption (character). Added with
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html).

## Value

A patchwork object (a ggplot-like object) that can be printed or saved.

## Details

Typical usage:

      decomp_plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
      grid_spr_decomp(decomp_plots, order = c("contrib","removed","survivorship"), ncol = 1)
