# tier4tools <img src="man/figures/logo.png" align="right" alt="" width="200" />
### Tier 4 Spawning Potential Ratio (SPR) reference points in Alaska groundfish stock assessments

`tier4tool` is an R package developed to support assessment scientists develop Tier 4 alternatives for groundfish aligned with North Pacific Fishery Management Council (NPFMC) guidance. The package provides a transparent, reproducible workflow for computing SPR-based reference points when full age-structured models are not available or not appropriate.

This package provides:
-  consistent implementation of SPR analyses for single-species and multispecies stocks,
-  report-ready visualizations and output,
-  diagnostic tools for understanding sources of uncertainty in SPR-based reference points.

## Installation and Development

`tier4tools` is under active development and is not yet on CRAN. You can install the development version from GitHuB:

```
# install.packages("remotes")
remotes::install_github("noaa-afsc/tier4tools")
```

Issues, suggestions, and contributions are welcome via GitHub.

## Basic Workflow

A typical Tier analysis using `tier4tools` follows four main steps:

1.  Define biological and fishery inputs using `spr_input()`. Visualize and validate inputs using `plot_spr_inputs()`.

2.  Compute per-recruit objects and SPR-based reference points (e.g., $F_{40\%}$, $F_{35\%}$) using `run_spr()`. Visualize results with `plot_spr_curves()`. 

3.  Show age decomposition and plus group diagnostics for spawning biomass per recruit (SBPR) contributions, removals, and survivorship using `plot_spr_decomp()`.

4.  Run structured sensitivity analyses on key SPR input parameters (e.g., maturity, fishery selectivity, natural mortality) using `spr_sensitivity()`. Visualize results using `plot_sensitivity_heatmap()`. For multispecies applications, sensitivity analyses for recruitment allocation by species are also available. These results are visualized using `plot_recprop_sensitivity()`.

## Vignettes

Detailed examples are available the package vignettes, including:

-  single-species Tier 4 workflows,
-  multispecies per-recruit analyses,
-  interpretation of age decomposition plots and sensitivity results

```
browseVignettes("tier4tools")
```

