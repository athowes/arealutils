#' Fit constant Small Area Estimation model using `rstan`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @template sf
#' @param nsim_warm Number of warmup samples, passed to `rstan`.
#' @param nsim_iter Number of samples, passed to `rstan`.
#' @param chains Number of chains, each of which gets `nsim_warm + nsim_iter` samples, passed to `rstan`.
#' @param cores Number of cores, passed to `rstan`, defaults to `parallel::detectCores()`
#' @examples
#' constant_stan(mw, cores = 2)
#' @export
constant_stan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  dat <- list(
    n_obs = n_obs,
    n_mis = n_mis,
    ii_obs = array(ii_obs),
    ii_mis = array(ii_mis),
    n = nrow(sf),
    y_obs = sf$y[ii_obs],
    m = sf$n_obs
  )
  
  fit <- rstan::sampling(
    stanmodels$constant,
    data = dat,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )
  
  return(fit)
}

#' Fit IID Small Area Estimation model using `rstan`.
#'
#' Random effects are independent and identically distributed.
#'
#' @inheritParams constant_stan
#' @examples
#' iid_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
iid_stan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)

  dat <- list(
    n_obs = n_obs,
    n_mis = n_mis,
    ii_obs = array(ii_obs),
    ii_mis = array(ii_mis),
    n = nrow(sf),
    y_obs = sf$y[ii_obs],
    m = sf$n_obs
  )

  fit <- rstan::sampling(
    stanmodels$iid,
    data = dat,
    warmup = nsim_warm,
    iter = nsim_iter,
    chains = chains,
    cores = cores
  )

  return(fit)
}

#' Fit Besag Small Area Estimation model using `rstan`.
#'
#' Random effects have an improper conditional autoregressive (ICAR)
#' distribution with (generalised) precision matrix produced using
#' the [`nb_to_precision`] function with input `nb`,
#' the neighbourhood structure of `sf`.
#'
#' @inheritParams constant_stan
#' @param method One of `"default"` or `"morris"`.
#' @examples
#' besag_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
besag_stan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), method = "default"){
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  scale <- get_scale(Q)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs)
  
  if(method == "default") {
    
    warning("Does not impose correct constraints!")
    
    Q_scaled <- scale_gmrf_precision(Q)$Q
    
    dat <- rlist::list.append(dat,
                              mu = rep(0, nrow(sf)),
                              Q = Q_scaled
    )
    
    fit <- rstan::sampling(
      stanmodels$mvn_precision,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  if(method == "morris") {
    
    warning("Assumes that the adjacency graph is connected!")
    
    g <- nb_to_graph(nb)
    
    dat <- rlist::list.append(dat,
                              n_edges = g$n_edges,
                              node1 = g$node1,
                              node2 = g$node2,
                              scaling_factor = scale
    )
    
    fit <- rstan::sampling(
      stanmodels$besag_morris,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}

#' Fit BYM2 Small Area Estimation model using `rstan`.
#'
#' @inheritParams constant_stan
#' @param method One of `"default"` or `"morris"`.
#' @examples
#' bym2_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
bym2_stan <- function(sf, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), method = "default"){
  
  warning("Doesn't take non-connectedness into account correctly!")
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  nb <- sf_to_nb(sf)
  Q <- nb_to_precision(nb)
  scale <- get_scale(Q)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs)
  
  if(method == "default") {
    
    warning("Does not impose correct constraints!")
    
    Q_scaled <- scale_gmrf_precision(Q)$Q
    
    dat <- rlist::list.append(dat,
                              mu = rep(0, nrow(sf)),
                              Q = Q_scaled
    )
    
    fit <- rstan::sampling(
      stanmodels$bym2_precision,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  if(method == "morris") {
    
    warning("Assumes that the adjacency graph is connected!")
    
    g <- nb_to_graph(nb)
    
    dat <- rlist::list.append(dat,
                              n_edges = g$n_edges,
                              node1 = g$node1,
                              node2 = g$node2,
                              scaling_factor = scale
    )
    
    fit <- rstan::sampling(
      stanmodels$bym2_morris,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}

#' Fit Centroid MVN Small Area Estimation model using `rstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyper-parameters
#' are fixed.
#'
#' @inheritParams constant_stan
#' @param bym2 Logical indicating if the centroid spatial random effects should be convoluted 
#' with unstructured IID noise, defaults to `FALSE`.
#' @inheritParams centroid_covariance
#' @examples
#' fck_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
fck_stan <- function(sf, bym2 = FALSE, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), kernel = matern, ...){
  
  cov <- centroid_covariance(sf, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs,
              Sigma = cov,
              mu = rep(0, nrow(sf)))
  
  if(bym2){
    fit <- rstan::sampling(
      stanmodels$bym2_covariance,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  else{
    fit <- rstan::sampling(
      stanmodels$mvn_covariance,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}

#' Fit Integrated MVN Small Area Estimation model using `rstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`].
#'
#' @inheritParams constant_stan
#' @param bym2 Logical indicating if the spatial random effects should be convoluted 
#' with unstructured IID noise, defaults to `FALSE`.
#' @inheritParams integrated_covariance
#' @examples
#' fik_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
fik_stan <- function(sf, bym2 = FALSE, L = 10, type = "hexagonal", nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores(), kernel = matern, ...){
  
  cov <- integrated_covariance(sf,  L = L, type = type, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs,
              Sigma = cov,
              mu = rep(0, nrow(sf)))
  
  if(bym2){
    fit <- rstan::sampling(
      stanmodels$bym2_centroid,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  else{
    fit <- rstan::sampling(
      stanmodels$centroid,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}

#' Fit Bayesian Centroid MVN Small Area Estimation model using `rstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_stan
#' @param bym2 Logical indicating if the spatial random effects should be convoluted 
#' with unstructured IID noise, defaults to `FALSE`.
#' @examples
#' ck_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ck_stan <- function(sf, bym2 = FALSE, nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
  
  D <- centroid_distance(sf)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(D)), plb = 0.01, pub = 0.01)
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs,
              mu = rep(0, nrow(sf)),
              a = param$a,
              b = param$b,
              D = D)
  
  if(bym2){
    fit <- rstan::sampling(
      stanmodels$bym2_centroid,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  else{
    fit <- rstan::sampling(
      stanmodels$centroid,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}

#' Fit Bayesian Integrated MVN Small Area Estimation model using `rstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`]. Kernel hyperparameters
#' are given a prior and learnt.
#'
#' @inheritParams constant_stan
#' @inheritParams integrated_covariance
#' @param bym2 Logical indicating if the spatial random effects should be convoluted 
#' with unstructured IID noise, defaults to `FALSE`.
#' @examples
#' ik_stan(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
#' @export
ik_stan <- function(sf, bym2 = FALSE, L = 10, type = "hexagonal", nsim_warm = 100, nsim_iter = 1000, chains = 4, cores = parallel::detectCores()){
  
  n <- nrow(sf)
  samples <- sf::st_sample(sf, type = type, exact = TRUE, size = rep(L, n))
  S <- sf::st_distance(samples, samples)
  
  # Parameters of the length-scale prior
  param <- invgamma_prior(lb = 0.1, ub = max(as.vector(S)), plb = 0.01, pub = 0.01)
  
  # Data structure for unequal number of points in each area
  sample_index <- sf::st_intersects(sf, samples)
  sample_lengths <- lengths(sample_index)
  start_index <- sapply(sample_index, function(x) x[1])
  
  ii_obs <- which(!is.na(sf$y))
  ii_mis <- which(is.na(sf$y))
  n_obs <- length(ii_obs)
  n_mis <- length(ii_mis)
  
  dat <- list(n_obs = n_obs,
              n_mis = n_mis,
              ii_obs = array(ii_obs),
              ii_mis = array(ii_mis),
              n = nrow(sf),
              y_obs = sf$y[ii_obs],
              m = sf$n_obs,
              mu = rep(0, nrow(sf)),
              a = param$a,
              b = param$b,
              sample_lengths = sample_lengths,
              total_samples = sum(sample_lengths),
              start_index = start_index,
              S = S)
  
  if(bym2){
    fit <- rstan::sampling(
      stanmodels$bym2_integrated,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  else{
    fit <- rstan::sampling(
      stanmodels$integrated,
      data = dat,
      warmup = nsim_warm,
      iter = nsim_iter,
      chains = chains,
      cores = cores
    )
  }
  
  return(fit)
}