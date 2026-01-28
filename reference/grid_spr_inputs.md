# Arrange SPR input plots into a report-ready grid

Combine list output from
[`plot_spr_inputs()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_inputs.md)
into a multi-panel figure using patchwork. Because patchwork is in
Suggests, this function will only run when it is installed.

## Usage

``` r
grid_spr_inputs(
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
  [`plot_spr_inputs()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_inputs.md).

- order:

  Character vector giving the plot names (in `plots`) in the order to
  arrange. Defaults to `names(plots)`.

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

      inp_plots <- plot_spr_inputs(inp, panels = c("mat_selex","wt","len","M"))
      grid_spr_inputs(inp_plots, order = c("mat_selex","wt","len","M"), ncol = 1)
