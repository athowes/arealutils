# arealutils

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

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

### Guidelines for adding a new `TMB` model

See [here](https://rdrr.io/github/mlysy/TMBtools/f/vignettes/TMBtools.Rmd).

1. Create `model.hpp` in `arealutils/src/TMB`
2. Run `TMBtools::export_models()`
3. Run `devtools::install()` to install the package locally

### Guidelines for adding a new `rstan` model

See [here](https://mc-stan.org/rstantools/articles/minimal-rstan-package.html).

1. Create `model.stan` in `arealutils/inst/stan`
2. Run `pkgbuild::compile_dll()` to perform a fake R CMD install
4. Run `devtools::install()` to install the package locally

## Models

| Model                                   | `TMB`   | `aghq`  | `tmbstan` | `R-INLA` | `rstan` |
|:----------------------------------------|:--------|:--------|:----------|:---------|:--------|
| Constant                                | &check; | &check; | &check;   | &check;  | -       | 
| Independent and identically distributed | &check; | &check; | &check;   | &check;  | -       | 
| Besag                                   | &check; | &check; | &check;   | &check;  | -       | 
| Besag-York-Mollié  2                    | -       | -       | -         | &check;  | -       | 
| Centroid kernel (fixed lengthscale)     | -       | -       | -         | &check;  | -       | 
| Integrated kernel (fixed lengthscale)   | -       | -       | -         | &check;  | -       | 
| Centroid kernel                         | -       | -       | -         | -        | -       | 
| Integrated kernel                       | -       | -       | -         | -        | -       | 

