#' Fit constant Small Area Estimation model using `tmbstan`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param nsim_warm Number of warmup samples, passed to `rstan`.
#' @param nsim_iter Number of samples, passed to `rstan`.
#' @param chains Number of chains, each of which gets `nsim_warm + nsim_iter` samples, passed to `rstan`.
#' @param cores Number of cores, passed to `rstan`, defaults to `parallel::detectCores()`.
#' @param ii The (zero-indexed) indices of the observations held-out.
#' @examples
#' constant_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
constant_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), ii = NULL){
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
iid_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), ii = NULL){
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
besag_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), ii = NULL){
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit BYM2 Small Area Estimation model using `tmbstan`.
#'
#' @inheritParams constant_tmbstan
#' @examples
#' bym2_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
bym2_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), ii = NULL){
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
                logit_phi = 0,
                log_sigma_u = 0)
  
  obj <- TMB::MakeADFun(
    data = c(model = "bym2", dat),
    parameters = param,
    random = c("beta_0", "u", "w"),
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

#' Fit Centroid MVN Small Area Estimation model using `tmbstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyper-parameters
#' are fixed.
#'
#' @inheritParams constant_tmbstan
#' @inheritParams centroid_covariance
#' @examples
#' fck_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
fck_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), kernel = matern, ii = NULL, ...){
  
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit Integrated MVN Small Area Estimation model using `tmbstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`].
#'
#' @inheritParams constant_tmbstan
#' @inheritParams integrated_covariance
#' @examples
#' fik_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
fik_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), L = 10, type = "random", kernel = matern, ii = NULL, ...){
  
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit Bayesian Centroid MVN Small Area Estimation model using `tmbstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_tmbstan
#' @examples
#' ck_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ck_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), ii = NULL){
  D <- centroid_distance(sf)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(D)), plb = 0.01, pub = 0.01)
  
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit Bayesian Integrated MVN Small Area Estimation model using `tmbstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_tmbstan
#' @inheritParams integrated_covariance
#' @examples
#' ik_tmbstan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ik_tmbstan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), L = 10, type = "random", ii = NULL, ...){
  n <- nrow(sf)
  samples <- sf::st_sample(sf, type = type, exact = TRUE, size = rep(L, n))
  S <- sf::st_distance(samples, samples)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(S)), plb = 0.01, pub = 0.01)
  
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
  
  fit <- tmbstan::tmbstan(
    obj = obj,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}