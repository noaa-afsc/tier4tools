# tier4tools

### Tier 4 Spawning Potential Ratio (SPR) reference points in Alaska groundfish stock assessments

`tier4tools` is an R package developed to support stock assessment
scientists develop North Pacific Fishery Management Council (NPFMC) Tier
4 alternatives for Alaskan groundfish stocks. The package provides a
transparent, reproducible workflow for computing SPR-based reference
points when full age-structured models are not available or not
appropriate.

This package provides:

- consistent implementation of SPR analyses for single-species and
  multispecies stocks,

- report-ready visualizations and output,

- diagnostic tools for understanding sources of uncertainty in SPR-based
  reference points.

## Installation and Development

`tier4tools` is under active development and is not yet on CRAN. You can
install the development version from GitHub:

    devtools::install_github("noaa-afsc/tier4tools", build_vignettes = TRUE)

Issues, suggestions, and contributions are welcome via GitHub.

## Basic Workflow

A typical Tier 4 analysis using `tier4tools` follows four main steps:

1.  Define biological and fishery inputs using
    [`spr_input()`](https://noaa-afsc.github.io/tier4tools/reference/spr_input.md).
    Visualize and validate inputs using
    [`plot_spr_inputs()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_inputs.md).

2.  Compute per-recruit objects and SPR-based reference points (e.g.,
    F40%, F35%) using
    [`run_spr()`](https://noaa-afsc.github.io/tier4tools/reference/run_spr.md).
    Visualize results with
    [`plot_spr_curves()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_curves.md).

3.  Show age decomposition and plus group diagnostics for spawning
    biomass per recruit (SBPR) contributions, removals, and survivorship
    using
    [`plot_spr_decomp()`](https://noaa-afsc.github.io/tier4tools/reference/plot_spr_decomp.md).

4.  Run structured sensitivity analyses on key SPR input parameters
    (e.g., maturity, fishery selectivity, natural mortality) using
    [`spr_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/spr_sensitivity.md).
    Visualize results using
    [`plot_sensitivity_heatmap()`](https://noaa-afsc.github.io/tier4tools/reference/plot_sensitivity_heatmap.md).
    For multispecies applications, sensitivity analyses for recruitment
    allocation by species are also available. These results can be
    visualized using
    [`plot_recprop_sensitivity()`](https://noaa-afsc.github.io/tier4tools/reference/plot_recprop_sensitivity.md).

## Vignettes

Detailed examples are available in the package vignettes, including:

- single-species Tier 4 workflows,

- multispecies per-recruit analyses,

- interpretation of age decomposition plots and sensitivity results

&nbsp;

    browseVignettes("tier4tools")
