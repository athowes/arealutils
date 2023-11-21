#' Compute length-scale such that averagely distant points have correlation `p`.
#'
#' See Best (1999) "Bayesian models for spatially correlated disease and exposure data".
#' When the average distance between points is below `1`m returns value one.
#' 
#' @param D Matrix of distances between points.
#' @param kernel A kernel function, defaults to `matern`.
#' @param p A percentage, defaults to 0.01.
#' @param ... Additional arguments to `kernel`.
#' @export
best_average <- function(D, kernel = matern, p = 0.01, ...) {
  # as.numeric to avoid units issues
  m <- as.numeric(mean(D))
  
  # When points are all the same (such that there is less than a metre between them) return l = 1
  if(m < 1) return(1)
  
  l_opt <- stats::uniroot(
    f = function(l) kernel(m, l, ...) - p, 
    lower = 1,
    upper = as.numeric(max(D)) # Length-scale unlikely to be greater than maximum distance
  )
  return(l_opt$root)
}

#' Find parameters of the inverse Gamma distribution such that `plb` of the probability mass is below `lb`
#' and `pub` of the probability mass is above `ub`.
#'
#' @param lb Lower bound.
#' @param ub Upper bound.
#' @param plb Proportion of the mass below `lb`, defaults to 0.01.
#' @param pub Proportion of the mass above `lb`, defaults to 0.01.
#' @source From \href{https://github.com/paul-buerkner/brms/commit/524e738aeeb82e49c5338839c3e10113763b6de1}{brms} 
#' following \href{Michael Betancourt}{https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html}.
#' @export
invgamma_prior <- function(lb, ub, plb = 0.01, pub = 0.01) {
  .opt_fun <- function(x, lb, ub) {
    # On the log-scale to ensure positive
    x <- exp(x)
    y1 <- invgamma::pinvgamma(lb, x[1], x[2], log.p = TRUE)
    y2 <- invgamma::pinvgamma(ub, x[1], x[2], lower.tail = FALSE, log.p = TRUE)
    c(y1 - log(plb), y2 - log(pub))
  }

  opt_res <- nleqslv::nleqslv(
    c(0, 0), .opt_fun, lb = lb, ub = ub,
    control = list(allowSingular = TRUE)
  )
  
  res <- exp(opt_res$x)
  return(list(a = res[1], b = res[2]))
}
