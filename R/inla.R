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

#' Fit IID Small Area Estimation model using `R-INLA`.
#'
#' Random effects are independent and identically distributed.
#'
#' @inheritParams constant_inla
#' @examples
#' iid_inla(mw, cores = 2)
#' @export
iid_inla <- function(sf, verbose = FALSE, cores = parallel::detectCores()){
  m <- NULL
  
  dat <- list(
    id = 1:nrow(sf),
    y = sf$y,
    m = sf$n_obs
  )
  
  # sigma ~ N(0. 2.5^2); initial in terms of log(tau) so 0 corresponds to tau = 1
  tau_prior <- list(
    prec = list(prior = "logtnormal", param = c(0, 1/2.5^2), initial = 0, fixed = FALSE)
  )
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  formula <- y ~ 1 + f(id, model = "iid", hyper = tau_prior)
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    verbose = verbose,
    num.threads = cores
  )
  
  return(fit)
}

#' Fit Besag Small Area Estimation model using `R-INLA`.
#'
#' Random effects have an improper conditional autoregressive (ICAR)
#' distribution. This is implemented by `R-INLA` using the
#' option `model = "besag"`.
#'
#' @inheritParams constant_inla
#' @examples
#' besag_inla(mw, cores = 2)
#' @export
besag_inla <- function(sf, verbose = FALSE, cores = parallel::detectCores()){
  
  m <- NULL
  
  nb <- sf_to_nb(sf)
  
  dat <- list(id = 1:nrow(sf),
              y = sf$y,
              m = sf$n_obs)
  
  # sigma ~ N(0. 2.5^2); initial in terms of log(tau) so 0 corresponds to tau = 1
  tau_prior <- list(prec = list(prior = "logtnormal", param = c(0, 1/2.5^2),
                                initial = 0, fixed = FALSE))
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  # constr = TRUE is a sum-to-zero constraint else +/- constant to all leaves density unchanged
  formula <- y ~ 1 + f(id,
                       model = "besag",
                       graph = nb,
                       scale.model = TRUE,
                       constr = TRUE,
                       hyper = tau_prior)
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m,
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

#' Fit BYM2 Small Area Estimation model using `R-INLA`.
#'
#' @inheritParams constant_inla
#' @examples
#' bym2_inla(mw, cores = 2)
#' @export
bym2_inla <- function(sf, verbose = FALSE, cores = parallel::detectCores()){
  
  m <- NULL
  
  nb <- sf_to_nb(sf)
  
  dat <- list(id = 1:nrow(sf),
              y = sf$y,
              m = sf$n_obs)
  
  # sigma ~ N(0. 2.5^2); initial in terms of log(tau) so 0 corresponds to tau = 1
  tau_prior <- list(prec = list(prior = "logtnormal", param = c(0, 1/2.5^2), initial = 0, fixed = FALSE))
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  # phi ~ Beta(0.5, 0.5); inital in terms of logit(phi), so 0 corresponds to phi = 0.5
  phi_prior <- list(phi = list(prior = "logitbeta", param = c(0.5, 0.5), initial = 0, fixed = FALSE))
  
  formula <- y ~ 1 + f(id,
                       model = "bym2",
                       graph = nb,
                       scale.model = TRUE,
                       constr = TRUE,
                       hyper = list(tau_prior, phi_prior))
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m,
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

inla.rgeneric.wicar.model <- function(
    cmd = c("graph", "Q", "mu", "initial", 
            "log.norm.const", "log.prior", "quit"), 
    theta = NULL) {
  
  envir <- parent.env(environment())
  
  interpret.theta <- function() {
    return(list(prec = exp(theta[1L])))
  }
  
  graph <- function() {
    return(Q())
  }
  
  Q <- function() {
    p <- interpret.theta()
    return(p$prec * R) # This stays sparse
  }
  
  mu <- function() {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
  }
  
  log.prior <- function() {
    # Havard's Gamma precision prior
    p <- interpret.theta()
    val <- stats::dgamma(p$prec, shape = 1, rate = 1, log = TRUE) + theta[1L]
    return(val)
  }
  
  initial <- function() {
    return(1)
  }
  
  quit <- function() {
    return(invisible())
  }
  
  if(is.null(theta)) {
    theta <- initial()
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}

#' Fit Weighted ICAR Small Area Estimation model using `R-INLA`.
#'
#' @param R A sparse precision matrix of class `"sparseMatrix"`
#' @inheritParams constant_inla
#' @export
wicar_inla <- function(sf, R, verbose = FALSE, cores = parallel::detectCores()){
  
  m <- NULL
  
  dat <- list(id = 1:nrow(sf),
              y = sf$y,
              m = sf$n_obs)
  
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  
  wicar <- INLA::inla.rgeneric.define(inla.rgeneric.wicar.model, R = R)
  
  # Missing scale.model, constr, hyper, ...
  formula <- y ~ 1 + f(id, model = wicar)
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = TRUE, waic = TRUE,
                           cpo = TRUE, config = TRUE),
    verbose = verbose,
    num.threads = cores
  )

  return(fit)
}

#' Fit Centroid MVN Small Area Estimation model using `R-INLA`.
#'
#' Random effects have a multivariate Gaussian distribution with covariance
#' matrix calculated using [`centroid_covariance`].
#'
#' @inheritParams constant_inla
#' @inheritParams centroid_covariance
#' @examples
#' fck_inla(mw, cores = 2)
#' @export
fck_inla <- function(sf, verbose = FALSE, kernel = matern, cores = parallel::detectCores(), ...){
  
  m <- NULL
  
  cov <- centroid_covariance(sf, kernel, ...)
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
  
  fit <- INLA::inla(
    formula,
    family = "xbinomial",
    control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior,
    data = dat,
    Ntrials = m,
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
  
  fit <- INLA::inla(
    formula,
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
    inla.mode = "experimental"
  )
  
  return(fit)
}