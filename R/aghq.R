#' Fit constant Small Area Estimation model using `aghq`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param k Number of quadrature points per hyperparameter dimension.
#' @param ii The (zero-indexed) indices of the observations held-out. 
#' @examples
#' constant_aghq(mw)
#' @export
constant_aghq <- function(sf, k = 3, ii = NULL){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs,
    left_out = !is.null(ii),
    ii = ifelse(!is.null(ii), ii, 0)
  )
  
  param <- list(
    beta_0 = 0
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "constant", dat),
    parameters = param,
    random = c("beta_0"),
    DLL = "arealutils_TMBExports"
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
#' iid_aghq(mw)
#' @export
iid_aghq <- function(sf, k = 3, ii = NULL){
  dat <- list(
    n = nrow(sf),
    y = sf$y,
    m = sf$n_obs,
    left_out = !is.null(ii),
    ii = ifelse(!is.null(ii), ii, 0)
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
#' besag_aghq(mw)
#' @export
besag_aghq <- function(sf, k = 3, ii = NULL){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
              Q = Q,
              Qrank = as.integer(Matrix::rankMatrix(Q)))
  
  param <- list(
    beta_0 = 0,
    u = rep(0, nrow(sf)),
    log_sigma_u = 0
  )
  
  obj <- TMB::MakeADFun(
    data = c(model = "besag", dat),
    parameters = param,
    random = c("beta_0", "u"),
    DLL = "arealutils_TMBExports"
  )
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit BYM2 Small Area Estimation model using `aghq`.
#'
#' @inheritParams constant_aghq
#' @examples
#' bym2_aghq(mw)
#' @export
bym2_aghq <- function(sf, k = 3, ii = NULL){
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  Q <- as(Q, "dgTMatrix")
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit Centroid MVN Small Area Estimation model using `aghq`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyper-parameters
#' are fixed.
#'
#' @inheritParams constant_aghq
#' @inheritParams centroid_covariance
#' @examples
#' fck_aghq(mw)
#' @export
fck_aghq <- function(sf, k = 3, kernel = matern, ii = NULL, ...){
  
  cov <- centroid_covariance(sf, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit Integrated MVN Small Area Estimation model using `aghq`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`].
#'
#' @inheritParams constant_aghq
#' @inheritParams integrated_covariance
#' @examples
#' fik_tmb(mw)
#' @export
fik_aghq <- function(sf, k = 3, L = 10, type = "hexagonal", kernel = matern, ii = NULL, ...){
  
  cov <- integrated_covariance(sf,  L = L, type = type, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit Bayesian Centroid MVN Small Area Estimation model using `aghq`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_aghq
#' @examples
#' ck_aghq(mw)
#' @export
ck_aghq <- function(sf, k = 3, ii = NULL){
  D <- centroid_distance(sf)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(D)), plb = 0.01, pub = 0.01)
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}

#' Fit Bayesian Integrated MVN Small Area Estimation model using `aghq`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_aghq
#' @inheritParams integrated_covariance
#' @examples
#' ik_aghq(mw)
#' @export
ik_aghq <- function(sf, k = 3, L = 10, type = "hexagonal", ii = NULL, ...){
  n <- nrow(sf)
  samples <- sf::st_sample(sf, type = type, exact = TRUE, size = rep(L, n))
  S <- sf::st_distance(samples, samples)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(S)), plb = 0.01, pub = 0.01)
  
  # Data structure for unequal number of points in each area
  sample_index <- sf::st_intersects(sf, samples)
  sample_lengths <- lengths(sample_index)
  start_index <- sapply(sample_index, function(x) x[1])
  
  dat <- list(n = nrow(sf),
              y = sf$y,
              m = sf$n_obs,
              left_out = !is.null(ii),
              ii = ifelse(!is.null(ii), ii, 0),
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}