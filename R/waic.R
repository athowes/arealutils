#' Watanabe–Akaike (Widely-Applicable) information criterion generic.
#'
#' @param fit Fitted model.
#' @param ... Additional arguments passed to `waic`.
#' @export
waic <- function(fit, ...) {
  UseMethod("waic")
}

#' @rdname waic
#' @export
waic.inla <- function(fit, ...) {
  local_waic <- fit$waic$local.waic
  est <- sum(local_waic)
  se <- stats::sd(local_waic) * sqrt(length(local_waic))
  return(list(est = est, se = se))
}

#' @rdname waic
#' @export
waic.stanfit <- function(fit, ...) {
  log_lik <- loo::extract_log_lik(fit)
  waic <- loo::waic(log_lik)
  est <- waic$estimates["waic", "Estimate"]
  se <- waic$estimates["waic", "SE"]
  return(list(est = est, se = se)) 
}

#' @rdname waic
#' @export
waic.sdreport <- function(fit, ...) {
  log_lik <- as.matrix(fit$value[which(names(fit$value) == "log_lik")])
  waic <- loo::waic(log_lik)
  est <- waic$estimates["waic", "Estimate"]
  se <- waic$estimates["waic", "SE"]
  return(list(est = est, se = se))
}

#' #' @rdname waic
#' #' @export
#' waic.aghq <- function(fit, ...) {
#'   log_lik <- NA
#'   waic <- loo::waic(log_lik)
#'   est <- waic$estimates["waic", "Estimate"]
#'   se <- waic$estimates["waic", "SE"]
#'   return(list(est = est, se = se)) 
#' }