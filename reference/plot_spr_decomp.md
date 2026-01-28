# Plot SPR diagnostics, SBPR age decompositions

Creates report-ready diagnostic plots of age-specific decompositions of
spawning biomass per recruit (SBPR) under unfished conditions (F = 0)
and at a reference fishing mortality (F = Ftarget) chosen by
[`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
(default = F40).

## Usage

``` r
plot_spr_decomp(
  spr_out,
  species = NULL,
  drop_plus_from_plot = TRUE,
  include_plus_caption = TRUE,
  which_panels = c("contrib", "removed", "survivorship")
)
```

## Arguments

- spr_out:

  Output list from
  [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
  with `diagnostics = TRUE`.

- species:

  Character vector of species to include. Defaults to all species in
  `spr_out`.

- drop_plus_from_plot:

  Logical. If TRUE (default), drops the terminal age from the plotted
  panels (only when plus-group diagnostics are available).

- include_plus_caption:

  Logical. If TRUE (default), adds the plus-group share caption to
  `plots$contrib` when available.

- which_panels:

  Character vector, any of
  `c("contrib", "removed", "removed_prop", "survivorship")`.

## Value

A named list of ggplot objects. A compact summary table is attached as
`attr(out, "summary")`.

## Details

Returned object is a named list of ggplot objects:

- `plots$contrib`: SBPR contribution by age, comparing F0 vs Ftarget.

- `plots$removed`: SBPR removed by fishing at each age (F0 - Ftarget).

- `plots$removed_prop`: (not shown by default) SBPR removed by fishing
  at each age, normalized to the species' total unfished SBPR (F0).
  Useful for multispecies applications to compare fractional losses
  among species when absolute SBPR magnitudes differ.

- `plots$survivorship`: survivorship by age (per recruit), F0 vs
  Ftarget.

Plus-group treatment can strongly influence SBPR, especially for
long-lived species. When `drop_plus_from_plot = TRUE` (default), the
terminal age is excluded from the plots to avoid dominance by the plus
group, while the fraction of total SBPR attributable to the terminal age
is retained in `spr_out$diagnostics$plus_share` and can be appended as a
plot caption.

**Summary table.** The function also attaches a compact per-species
summary table (unfished SBPR, SBPR at Ftarget, implied SPR at Ftarget,
and plus-group shares when available) as `attr(plots, "summary")`.

## Examples

``` r
if (FALSE) { # \dontrun{
spr_out <- run_spr(inp, diagnostics = TRUE,
                  multispecies_constraint = "all_species")

# Default panels
plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
plots$contrib
plots$removed
plots$survivorship

# Summary table for optional reporting
attr(plots, "summary")
} # }
```
