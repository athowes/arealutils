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

  fit <- INLA::inla(formula,
                    family = "xbinomial",
                    control.family = list(control.link = list(model = "logit")),
                    control.fixed = beta_prior,
                    data = dat,
                    Ntrials = m,
                    control.predictor = list(compute = TRUE, link = 1),
                    control.compute = list(dic = TRUE, waic = TRUE,
                                           cpo = TRUE, config = TRUE),
                    verbose = verbose,
                    num.threads = cores)
  
  return(fit)
}