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

#' Fit constant Small Area Estimation model using `R-INLA`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @template sf
#' @param verbose Should `R-INLA` run in mode `verbose = TRUE`.
#' @param cores Number of cores, passed to `R-INLA`, defaults to `parallel::detectCores()`.
#' @examples
#' fit <- constant_inla(mw, cores = 2)
#' mean <- plogis(fit$summary.fixed$mean)
#' @export
constant_inla <- function(sf, verbose = FALSE, cores = parallel::detectCores()){
  m <- NULL
  
  dat <- list(
    id = 1:nrow(sf),
    y = sf$y,
    m = sf$n_obs
  )
  
  formula <- y ~ 1
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m, # Picks up the correct column in the dataframe dat
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(
      dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE, 
      return.marginals.predictor = TRUE
    ),
    verbose = verbose,
    num.threads = cores,
    inla.mode = "experimental"
  )
  
  return(fit)
}

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

#' Fit constant Small Area Estimation model using `aghq`.
#'
#' Simply fits a constant (the mean). This is useful as a benchmark for other models.
#'
#' @param sf A simple features object with some geometry.
#' @param its Number of iterations in outer loop optimisation, passed to `nlminb`.
#' @examples
#' constant_tmb(mw, its = 100)
#' @export
constant_aghq <- function(sf, its = 1000){
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
  
  quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
  
  return(quad)
}