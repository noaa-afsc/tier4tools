# Plot multispecies sensitivity to recruitment proportions

Creates report-ready plots summarizing how multispecies (complex-level)
SPR reference points change across recruitment proportion assumptions
for two-species complexes. The function is designed to work directly
with the output of
[`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md)
run under two settings: (1) an unconstrained run where complex reference
points are derived from the combined SPR curve
(`multispecies_constraint = "none"`), and (2) a constrained run where
all species must meet the target SPR
(`multispecies_constraint = "all_species"`).

## Usage

``` r
plot_recprop_sensitivity(
  sens_uncon,
  sens_con,
  targets = c(0.4, 0.35),
  prop_col = "prop_value",
  show_delta = FALSE,
  include_species = FALSE,
  include_yield_share = FALSE,
  title = "Sensitivity of multispecies reference points to recruitment proportions",
  xlab = NULL
)
```

## Arguments

- sens_uncon:

  A data frame returned by
  [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md)
  with `multispecies_constraint = "none"`.

- sens_con:

  A data frame returned by
  [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md)
  with `multispecies_constraint = "all_species"`.

- targets:

  Numeric vector of SPR targets to plot (values in (0, 1)), for example
  `c(0.40, 0.35)`.

- prop_col:

  Character. Name of the recruitment proportion column to use on the
  x-axis. Defaults to `"prop_value"` (recommended), but can be set to
  `"prop"` if desired.

- show_delta:

  Logical. If FALSE (default), plots absolute reference point F. If
  TRUE, plots deviations from the unconstrained reference point
  (\\\Delta F\\).

- include_species:

  Logical. If TRUE, attempts to return a limiting-species plot as a
  second panel.

- include_yield_share:

  Logical. If TRUE, attempts to return an implied yield-share plot as a
  third panel.

- title:

  Character. Main plot title.

- xlab:

  Character. Optional x-axis label. If `NULL`, the function uses the
  species names stored in `attr(sens_uncon, "species")` when available.

## Value

If `include_species = FALSE` and `include_yield_share = FALSE`, returns
a single ggplot object (the main plot).

Otherwise returns a named list with:

- `main`:

  The main reference-point sensitivity plot (ggplot).

- `limiting_species`:

  A limiting-species panel (ggplot) or `NULL`.

- `yield_share`:

  An implied yield-share panel (ggplot) or `NULL`.

## Details

**Default behavior.** The main plot is a single panel where:

- color = SPR target (for example, 0.40 vs 0.35), and

- linetype = the reference point method (unconstrained vs constrained).

**Delta option.** When `show_delta = TRUE`, the unconstrained values are
plotted as a 0 baseline and the constrained values shows the
constraint-induced change at each recruitment proportion and target:
\$\$\Delta F = F\_{\mathrm{method}} - F\_{\mathrm{unconstrained}}.\$\$

**Optional diagnostics.**

- When `include_species = TRUE`, the function also attempts to return a
  limiting-species panel based on species-specific reference point
  columns produced by
  [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md)
  (for example `F40_RE`, `F40_BS`). The limiting species is defined as
  the species with the smallest species-specific reference point at each
  recruitment proportion and target.

- When `include_yield_share = TRUE`, the function attempts to return an
  implied yield-share panel using yield-share columns added by
  [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md),
  for example `share40_RE`, `share40_BS`, and constrained equivalents
  `share40_RE_constrained`, `share40_BS_constrained`. Yield shares are
  interpreted as recruitment-weighted YPR contributions at the reference
  point F.

If the required columns are not present, the corresponding optional
plots are returned as `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
sens_uncon <- spr_sensitivity(inp, prop = seq(0.2, 0.8, 0.05),
                             multispecies_constraint = "none")
sens_con   <- spr_sensitivity(inp, prop = seq(0.2, 0.8, 0.05),
                             multispecies_constraint = "all_species")

# Main plot only
plot_recprop_sensitivity(sens_uncon, sens_con)

# Delta plot
plot_recprop_sensitivity(sens_uncon, sens_con, show_delta = TRUE)

# Return diagnostic panels (if columns are available)
plots <- plot_recprop_sensitivity(
  sens_uncon, sens_con,
  include_species = TRUE,
  include_yield_share = TRUE
)
plots$main
plots$limiting_species
plots$yield_share
} # }
```
