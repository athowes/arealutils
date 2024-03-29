#' Fit constant Small Area Estimation model using `TMB`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param its Number of iterations in outer loop optimisation, passed to `nlminb`.
#' @param ii The (zero-indexed) indices of the observations held-out.
#' @examples
#' constant_tmb(mw, its = 100)
#' @export
constant_tmb <- function(sf, its = 1000, ii = NULL){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs,
    left_out = !is.null(ii),
    ii = if(!is.null(ii)) ii else 0
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
#' @inheritParams constant_tmb
#' @examples
#' iid_tmb(mw, its = 100)
#' @export
iid_tmb <- function(sf, its = 1000, ii = NULL){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs,
    left_out = !is.null(ii),
    ii = if(!is.null(ii)) ii else 0
  )
  
  param <- list(
    beta_0 = 0,
    u = rep(0, nrow(sf)),
    log_sigma_u = 0
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "iid", dat),
    parameters = param,
    random = c("beta_0", "u"),
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
#' distribution with (generalised) precision matrix produced using
#' the [`nb_to_precision`] function with input `nb`,
#' the neighbourhood structure of `sf`.
#'
#' @inheritParams constant_tmb
#' @examples
#' besag_tmb(mw, its = 100)
#' @export
besag_tmb <- function(sf, its = 1000, ii = NULL){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              Q = Q,
              Qrank = as.integer(Matrix::rankMatrix(Q)))
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                log_sigma_u = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "besag", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit BYM2 Small Area Estimation model using `TMB`.
#'
#' @inheritParams constant_tmb
#' @examples
#' bym2_tmb(mw, its = 100)
#' @export
bym2_tmb <- function(sf, its = 1000, ii = NULL){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              Q = Q)
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                w = rep(0, dat$n),
                log_sigma_u = 0,
                logit_phi = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "bym2", dat),
    parameters = param,
    random = c("beta_0", "u", "w"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit Centroid MVN Small Area Estimation model using `TMB`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyper-parameters
#' are fixed.
#'
#' @inheritParams constant_tmb
#' @inheritParams centroid_covariance
#' @examples
#' fck_tmb(mw, its = 100)
#' @export
fck_tmb <- function(sf, its = 1000, kernel = matern, ii = NULL, ...){
  
  cov <- centroid_covariance(sf, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              Sigma = cov)
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                log_sigma_u = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "mvn_covariance", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit Integrated MVN Small Area Estimation model using `TMB`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`].
#'
#' @inheritParams constant_tmb
#' @inheritParams integrated_covariance
#' @examples
#' fik_tmb(mw, its = 100)
#' @export
fik_tmb <- function(sf, its = 1000, L = 10, type = "hexagonal", kernel = matern, ii = NULL, ...){
  
  cov <- integrated_covariance(sf,  L = L, type = type, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              Sigma = cov)
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                log_sigma_u = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "mvn_covariance", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit Bayesian Centroid MVN Small Area Estimation model using `TMB`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_tmb
#' @examples
#' ck_tmb(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ck_tmb <- function(sf, its = 1000, ii = NULL){
  D <- centroid_distance(sf)
  
  # Parameters of the length-scale prior
  if(inherits(D, "units")) {
    if(units(D)$numerator == "m") {
      D <- units::set_units(D, "km")
    }
  }
  
  D_nonzero <- as.vector(D)[as.vector(D) > 0]
  lb <- quantile(D_nonzero, 0.05)
  ub <- quantile(D_nonzero, 0.95)
  
  param <- invgamma_prior(lb = lb, ub = ub, plb = 0.05, pub = 0.05)
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              a = param$a,
              b = param$b,
              D = D)
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                log_sigma_u = 0,
                log_l = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "centroid", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}

#' Fit Bayesian Integrated MVN Small Area Estimation model using `TMB`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_tmb
#' @inheritParams integrated_covariance
#' @examples
#' ik_tmb(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ik_tmb <- function(sf, its = 1000, L = 10, type = "hexagonal", ii = NULL, ...){
  n <- nrow(sf)
  samples <- sf::st_sample(sf, type = type, exact = TRUE, size = rep(L, n))
  S <- sf::st_distance(samples, samples)
  
  # Parameters of the length-scale prior
  if(inherits(S, "units")) {
    if(units(S)$numerator == "m") {
      S <- units::set_units(S, "km")
    }
  }
  
  S_nonzero <- as.vector(S)[as.vector(S) > 0]
  lb <- quantile(S_nonzero, 0.05)
  ub <- quantile(S_nonzero, 0.95)
  
  param <- invgamma_prior(lb = lb, ub = ub, plb = 0.05, pub = 0.05)
  
  # Data structure for unequal number of points in each area
  sample_index <- sf::st_intersects(sf, samples)
  sample_lengths <- lengths(sample_index)
  
  # Take one away due to zero indexing in C++
  start_index <- sapply(sample_index, function(x) x[1] - 1)
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = if(!is.null(ii)) ii else 0,
              a = param$a,
              b = param$b,
              sample_lengths = sample_lengths,
              total_samples = sum(sample_lengths),
              start_index = start_index,
              S = S)
  
  param <- list(beta_0 = 0,
                u = rep(0, dat$n),
                log_sigma_u = 0,
                log_l = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "integrated", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  opt <- nlminb(start = obj$par,
                objective = obj$fn,
                gradient = obj$gr,
                control = list(iter.max = its, trace = 0))
  
  sd_out <- TMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE)
  
  return(sd_out)
}