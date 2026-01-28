# Plot SPR–F curves with reference points

Visualize spawning potential ratio (SPR) as a function of fishing
mortality (F) from a
[`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
output, with report-ready annotations of SPR-based fishing mortality
reference points (for example, F40% and F35%).

## Usage

``` r
plot_spr_curves(spr_out, which = c("total", "species"), digits = 3)
```

## Arguments

- spr_out:

  Output list from
  [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md).

- which:

  Character string. For multispecies outputs only:

  `"total"`

  :   Plot combined SPR(F) for the complex (black), with species
      overlays (colored).

  `"species"`

  :   Plot species-specific SPR(F) curves (colored).

- digits:

  Integer number of digits to show in reported F values.

## Value

A ggplot object.

## Details

**What is plotted.** The function plots `SPR(F) = SBPR(F) / SBPR(0)`
across the fishing mortality grid `spr_out$F_vec`. Horizontal reference
lines are drawn at each target SPR level in `spr_out$spr_targets` and
vertical reference lines are drawn at each F-based reference points.

**Single species.** If `spr_out` contains only one species, the function
always produces the single-species plot and ignores `which`.

**Multispecies complexes.**

- `which = "species"` shows one SPR(F) curve per species (colored) and
  annotates each species' own reference points.

- `which = "total"` shows the combined complex SPR(F) curve (black) with
  optional species overlays (colored). If constrained complex reference
  points exist (i.e., `spr_out$F_spr_total_constrained` contains finite
  values), the plot also shows both unconstrained and constrained
  complex reference points, and highlights the limiting species at the
  constrained reference points.

**Multispecies constraints.** If
[`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
was called with `multispecies_constraint = "all_species"`, constrained
combined reference points may be available in
`spr_out$F_spr_total_constrained`. When present, constrained and
unconstrained reference points are shown together and distinguished by
`Method` (shape and linetype).

## Examples

``` r
if (FALSE) { # \dontrun{
spr_out1 <- run_spr(inp_single, diagnostics = FALSE)
plot_spr_curves(spr_out1)

spr_out2 <- run_spr(inp_multi, diagnostics = FALSE, multispecies_constraint = "all_species")
plot_spr_curves(spr_out2, which = "total")
plot_spr_curves(spr_out2, which = "species")
} # }
```
