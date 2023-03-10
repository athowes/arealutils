#' Fit Bayesian Centroid MVN Small Area Estimation model using `rstan`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`]. Unlike `m5_stan`, this
#' version has a fully Bayesian treatment of kernel hyperparameters.
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
    fit <- rstan::sampling(stanmodels$bym2_centroid,
                           data = dat,
                           warmup = nsim_warm,
                           iter = nsim_iter,
                           chains = chains,
                           cores = cores)
  }
  else{
    fit <- rstan::sampling(stanmodels$centroid,
                           data = dat,
                           warmup = nsim_warm,
                           iter = nsim_iter,
                           chains = chains,
                           cores = cores)
  }

  return(fit)
}