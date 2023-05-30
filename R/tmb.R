#' Fit constant Small Area Estimation model using `TMB`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param its Number of iterations in outer loop optimisation, passed to `nlminb`.
#' @examples
#' constant_tmb(mw, its = 100)
#' @export
constant_tmb <- function(sf, its = 1000){
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
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit IID Small Area Estimation model using `TMB`.
#'
#' Random effects are independent and identically distributed.
#'
#' @param sf A simple features object with some geometry.
#' @param its Number of iterations in outer loop optimisation, passed to `nlminb`.
#' @examples
#' iid_tmb(mw, its = 100)
#' @export
iid_tmb <- function(sf, its = 1000){
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
  
  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(iter.max = its, trace = 0)
  )
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit Besag Small Area Estimation model using `TMB`.
#'
#' Random effects have an improper conditional autoregressive (ICAR)
#' distribution.
#'
#' @param sf A simple features object with some geometry.
#' @param its Number of iterations in outer loop optimisation, passed to
#' \code{nlminb}.
#' @examples
#' besag_tmb(mw, its = 100)
besag_tmb <- function(sf, its = 1000){
  
  compile("tmb/model2.cpp")
  dyn.load(dynlib("tmb/model2"))
  
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              Q = Q)
  
  param <- list(beta_0 = 0,
                phi = rep(0, dat$n),
                l_sigma_phi = 0)
  
  obj <- MakeADFun(data = dat,
                   parameters = param,
                   DLL = "model2")
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- sdreport(obj,
                     par.fixed = opt$par,
                     getJointPrecision = TRUE)
  
  return(sd_out)
}