# arealutils

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

<!-- badges: end -->

## Installation

This R package is being developed in support of the analysis in [`beyond-borders`](https://github.com/athowes/beyond-borders).
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
| Constant                                | &check; | &check; | &check;   | &check;  | &rarr;  | 
| Independent and identically distributed | &check; | &check; | &check;   | &check;  | &rarr;  | 
| Besag                                   | &check; | &check; | &check;   | &check;  | &rarr;  |  
| Besag-York-Mollié  2                    | &check; | &check; | &check;   | &check;  | &rarr;  |  
| Centroid kernel (fixed lengthscale)     | &check; | &check; | &check;   | &check;  | &rarr;  |  
| Integrated kernel (fixed lengthscale)   | &check; | &check; | &check;   | &check;  | &rarr;  |  
| Centroid kernel                         | &check; | &check; | &check;   | -        | &rarr;  |  
| Integrated kernel                       | &check; | &check; | &check;   | -        | &rarr;  | 


* &check;: done
* &rarr;: done but temporarily stored elsewhere (`rstan` compilation times when building the package slow development)

## Common tasks

* To build the website, use `pkgdown::build_site()`
* To document functions, use `devtools::document()`

## Initialisation

The R code used to initialise this package (much of it perhaps not required!) was:

```r
package_name <- "arealutils"

TMBtools::tmb_create_package(
  path = package_name,
  example_code = TRUE,
  fields = list(
  `Authors@R` = 'person("Your", "Name", email = "valid@email.com", role = c("aut", "cre"))')
)

# R wrapper functions to TMB models
usethis::use_template(
  template = "norm_ADFun.R",
  package = "TMBtools",
  save_as = file.path("R", "norm_ADFun.R"),
  data = list(pkg = package_name)
)

usethis::use_template(
  template = "gamma_ADFun.R",
  package = "TMBtools",
  save_as = file.path("R", "gamma_ADFun.R"),
  data = list(pkg = package_name)
)

# testthat tests
usethis::use_testthat()

usethis::use_package(package = "numDeriv", type = "Suggests")

usethis::use_template(
  template = "test-norm_ADFun.R",
  package = "TMBtools",
  save_as = file.path("tests", "testthat", "test-norm_ADFun.R")
)

usethis::use_template(
  template = "test-gamma_ADFun.R",
  package = "TMBtools",
  save_as = file.path("tests", "testthat", "test-gamma_ADFun.R")
)

# update NAMESPACE and package documentation
pkgbuild::compile_dll() # need to compile src first
devtools::document()

rstantools::use_rstan(pkgdir = ".", license = TRUE, auto_config = TRUE)
rstantools::rstan_config(pkgdir = ".")
```
