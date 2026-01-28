# Sensitivity analysis for SPR reference points

Runs a grid-based sensitivity analysis for Tier 4 SPR calculations by
varying a user-defined set of influential inputs. Available options
include: scalar natural mortality (M), logistic maturity \\a\_{50}\\,
knife-edge selectivity \\a_c\\, and (for 2-species complexes)
recruitment proportions.

## Usage

``` r
spr_sensitivity(
  x_base,
  M = NULL,
  maturity_a50 = NULL,
  selex_ac = NULL,
  prop = NULL,
  F_vec = seq(0, 0.4, by = 0.002),
  spr_targets = c(0.4, 0.35),
  multispecies_constraint = c("none", "all_species"),
  use_constrained_for_yield = TRUE,
  digits = 3
)
```

## Arguments

- x_base:

  A
  [`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
  object used as the base scenario.

- M:

  Numeric vector of scalar M values to apply. If `NULL`, M is not
  varied.

- maturity_a50:

  Numeric vector of logistic maturity a50 values. If `NULL`, not varied.

- selex_ac:

  Numeric vector of knife-edge selectivity a_c values. If `NULL`, not
  varied.

- prop:

  Numeric vector of recruitment proportions for 2-species inputs.
  Interpreted as the proportion assigned to the first species in
  `names(x_base$species)`. Ignored for 1 species.

- F_vec:

  Fishing mortality grid passed to
  [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md).

- spr_targets:

  SPR targets passed to
  [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
  (default `c(0.40, 0.35)`).

- multispecies_constraint:

  Passed to
  [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
  (`"none"` or `"all_species"`).

- use_constrained_for_yield:

  Logical. If TRUE and constrained F40 exists, evaluate `YPR_total` at
  constrained F40; otherwise evaluate at unconstrained F40.

- digits:

  Digits used for formatting `scenario_label`.

## Value

A df with one row per scenario. Columns include:

- `scenario_id`:

  Unique integer id.

- `M`, `maturity_a50`, `selex_ac`, `prop`:

  Scenario inputs (NA means base).

- `F40_total`, `F35_total`:

  Combined reference points (unconstrained).

- `F40_total_constrained`, `F35_total_constrained`:

  Combined constrained reference points when
  `multispecies_constraint="all_species"` (otherwise NA).

- `YPR_at_F40`:

  Total yield per recruit evaluated at F40 (constrained if requested).

- `plus_share_<sp>_F0`:

  Plus-group share at F=0 for each species `<sp>`.

- `plus_share_<sp>_Ftarget`:

  Plus-group share at the diagnostic Ftarget for each species.

- `scenario_label`:

  Readable label concatenating scenario inputs.

## Details

This function modifies a base
[`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md)
object `x_base` by updating one or more schedules and re-running
[`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md)
for each scenario. The result is returned as a tidy data frame suitable
for sensitivity plots or tables.

**Important limitation:** All sensitivity sweep dimensions are applied
to *all species simultaneously*. For example, sweeping `maturity_a50`
sets the maturity ogive midpoint to the same value for every species in
the input, and sweeping `selex_ac` sets the same knife-edge \\a_c\\ for
every species. Species-specific sensitivity sweeps (varying a parameter
for one species while holding others fixed) are not currently
implemented.

**What can be varied:**

- `M`: scalar \\M\\ applied at all ages (replaces existing `M_at_age`).

- `maturity_a50`: maturity ogive midpoint \\a\_{50}\\ for a logistic
  maturity schedule.

- `selex_ac`: knife-edge selectivity age \\a_c\\.

- `prop`: recruitment proportion assigned to the first species
  (2-species inputs only).

**Outputs:** The returned df includes combined (total) and (optionally)
constrained SPR reference points (e.g., F40, F35), total YPR evaluated
at F40, and plus group diagnostics as separate columns for each species.

## Examples

``` r
if (FALSE) { # \dontrun{
out <- spr_sensitivity(
  x_base,
  M = seq(0.03, 0.06, by = 0.005),
  maturity_a50 = seq(10, 20, by = 2),
  selex_ac = c(8, 10, 12),
  prop = seq(0.3, 0.7, by = 0.1),
  F_vec = seq(0, 0.4, by = 0.002),
  multispecies_constraint = "all_species"
)

# Heatmap: F40_total across (M, a50)
plot_sensitivity_heatmap(out, x = "M", y = "maturity_a50", fill = "F40_total")
} # }
```
