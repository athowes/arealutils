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
#' constant_tmbstan(mw, its = 100)
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