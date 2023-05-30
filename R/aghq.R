#' Fit constant Small Area Estimation model using `aghq`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param k Number of quadrature points per hyperparameter dimension.
#' @param its Number of iterations in outer loop optimisation, passed to `nlminb`.
#' @examples
#' constant_aghq(mw, its = 100)
#' @export
constant_aghq <- function(sf, k = 3, its = 1000){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs
  )
  
  param <- list(
    beta_0 = 0
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "constant", dat),
    parameters = param,
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = its, trace = 0)
  )
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit IID Small Area Estimation model using `aghq`.
#'
#' Random effects are independent and identically distributed.
#'
#' @inheritParams constant_aghq
#' @examples
#' iid_aghq(mw, its = 100)
#' @export
iid_aghq <- function(sf, k = 3, its = 1000){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs
  )
  
  param <- list(
    beta_0 = 0,
    phi = rep(0, nrow(sf)),
    sigma_phi = 1
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "iid", dat),
    parameters = param,
    random = c("beta_0", "phi"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = its, trace = 0)
  )
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit Besag Small Area Estimation model using `aghq`.
#'
#' Random effects have an improper conditional autoregressive (ICAR)
#' distribution with (generalised) precision matrix produced using
#' the [`nb_to_precision`] function with input `nb`,
#' the neighbourhood structure of `sf`.
#'
#' @inheritParams constant_aghq
#' @examples
#' besag_aghq(mw, its = 100)
#' @export
besag_aghq <- function(sf, k = 3, its = 1000){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs
  )
  
  param <- list(
    beta_0 = 0,
    phi = rep(0, nrow(sf)),
    sigma_phi = 1
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "iid", dat),
    parameters = param,
    random = c("beta_0", "phi"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = its, trace = 0)
  )
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit BYM2 Small Area Estimation model using `aghq`.
#'
#' @inheritParams constant_aghq
#' @examples
#' bym2_aghq(mw, its = 100)
#' @export
bym2_aghq <- function(sf, k = 3, its = 1000){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              Q = Q)
  
  param <- list(beta_0 = 0,
                phi = rep(0, dat$n),
                u = rep(0, dat$n),
                logit_pi = 0,
                sigma_phi = 1)
  
  obj <- TMB::MakeADFun(
    data = c(model = "bym2", dat),
    parameters = param,
    random = c("beta_0", "phi", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}