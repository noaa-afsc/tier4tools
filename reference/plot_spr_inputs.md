# Plot SPR input assumptions

Visualize inputs supplied to spr_input() object for checking assumptions
and reporting (base vs sensitivity, or across species). By default,
makes:

- maturity and selectivity on the same panel

- weight-at-age Optional panels include length-at-age and M-at-age.

## Usage

``` r
plot_spr_inputs(
  x,
  species = NULL,
  compare = NULL,
  panels = c("mat_selex", "wt", "len", "M"),
  overlay_maturity_selex = TRUE,
  show_fleets = FALSE,
  return_data = FALSE
)
```

## Arguments

- x:

  An object returned by spr_input().

- species:

  Character vector of species names to plot. Defaults to all in
  x\$species.

- compare:

  Optional named list of additional spr_input objects (e.g.,
  sensitivities). Each element name is used as the scenario label in
  legends/facets. Example: compare = list(Base = inp_base, Sens1 =
  inp_sens)

- panels:

  Character vector specifying which panels to draw. Any of:
  c("mat_selex","wt","len","M")

- overlay_maturity_selex:

  Logical. If TRUE, maturity and selectivity are overlayed on a single
  plot (recommended).

- show_fleets:

  Logical. If TRUE, and if fleet-specific selectivity is available,
  overlays fleet selectivity curves and the effective selectivity curve.
  Default FALSE.

- return_data:

  Logical. If TRUE, returns a list with plots and the tidy data used.

## Value

By default, returns a named list of ggplot objects. If return_data=TRUE,
returns list(plots=..., data=...).
