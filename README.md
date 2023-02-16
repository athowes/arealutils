# arealutils

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R build
status](https://github.com/athowes/arealutils/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/athowes/arealutils/actions)

<!-- badges: end -->

## Installation

This R package is being developed in support of the analysis in [`areal-comparison`](https://github.com/athowes/areal-comparison).
For this reason, many of the functions are relatively unlikely to be of great general use.
This disclaimer given, the package can be installed via:

```r
devtools::install_github("athowes/arealutils")
```

## `TMB` and `rstan`

This package uses compiled models from both `TMB` and `rstan`.
To do this we make use of the very helpful `TMBtools` and `rstantools` to set-up the package.
