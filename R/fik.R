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
    fit <- rstan::sampling(stanmodels$bym2_covariance,
                           data = dat,
                           warmup = nsim_warm,
                           iter = nsim_iter,
                           chains = chains,
                           cores = cores)
  }
  else{
    fit <- rstan::sampling(stanmodels$mvn_covariance,
                           data = dat,
                           warmup = nsim_warm,
                           iter = nsim_iter,
                           chains = chains,
                           cores = cores)
  }
  
  return(fit)
}

#' Fit Sampling MVN Small Area Estimation model using `R-INLA`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`integrated_covariance`].
#'
#' @inheritParams constant_inla
#' @inheritParams integrated_covariance
#' @examples
#' fik_inla(mw, cores = 2)
#' @export
fik_inla <- function(sf, verbose = FALSE, L = 50, kernel = matern, cores = parallel::detectCores(), ...){
  
  m <- NULL
  
  cov <- integrated_covariance(sf, L, kernel, ...)
  cov <- cov / riebler_gv(cov) # Standardise so tau prior is right
  C <- Matrix::solve(cov) # Precision matrix
  
  dat <- list(id = 1:nrow(sf),
              y = sf$y,
              m = sf$n_obs)
  
  # sigma ~ N(0. 2.5^2); initial in terms of log(tau) so 0 corresponds to tau = 1
  tau_prior <- list(prec = list(prior = "logtnormal", param = c(0, 1/2.5^2),
                                initial = 0, fixed = FALSE))
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  # See inla.doc("generic0")
  formula <- y ~ 1 + f(id, model = "generic0", Cmatrix = C, hyper = tau_prior)
  
  fit <- INLA::inla(formula,
                    family = "xbinomial",
                    control.family = list(control.link = list(model = "logit")),
                    data = dat,
                    control.fixed = beta_prior,
                    Ntrials = m,
                    control.predictor = list(compute = TRUE, link = 1),
                    control.compute = list(
                      dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE, 
                      return.marginals.predictor = TRUE
                    ),
                    verbose = verbose,
                    num.threads = cores,
                    inla.mode = "experimental")
  
  return(fit)
}