#' Fit constant Small Area Estimation model using `tmbstan`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param nsim_warm Number of warmup samples, passed to `rstan`.
#' @param nsim_iter Number of samples, passed to `rstan`.
#' @param chains Number of chains, each of which gets `nsim_warm + nsim_iter` samples, passed to `rstan`.
#' @param cores Number of cores, passed to `rstan`, defaults to `parallel::detectCores()`
#' @examples
#' constant_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
constant_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit IID Small Area Estimation model using `tmbstan`.
#'
#' Random effects are independent and identically distributed.
#'
#' @inheritParams constant_tmbstan
#' @examples
#' constant_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
iid_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
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
    random = "phi",
    DLL = "arealutils_TMBExports"
  )
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit Besag Small Area Estimation model using `tmbstan`.
#'
#' Random effects have an improper conditional autoregressive (ICAR)
#' distribution with (generalised) precision matrix produced using
#' the [`nb_to_precision`] function with input `nb`,
#' the neighbourhood structure of `sf`.
#'
#' @inheritParams constant_tmbstan
#' @examples
#' besag_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
besag_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              Q = Q,
              Qrank = as.integer(Matrix::rankMatrix(Q)))
  
  param <- list(beta_0 = 0,
                phi = rep(0, dat$n),
                sigma_phi = 1)
  
  obj <- TMB::MakeADFun(
    data = c(model = "besag", dat),
    parameters = param,
    random = "phi",
    DLL = "arealutils_TMBExports"
  )
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}
