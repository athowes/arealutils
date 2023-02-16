#' Compute forecast assessments for model at a single held-out point.
#'
#' @template sf
#' @param fit Fitted model.
#' @param i The index (in `1:nrow(sf)`) of the held-out point.
#' @param S The number of Monte Carlo samples to draw from the approximate posterior.
#' @export
held_out_metrics <- function(fit, sf, i, S = 1000){
  y <- sf$y[[i]]
  n_obs <- sf$n_obs[[i]]
  s <- sample_marginal(fit, i, n_obs, S)
  error_samples <- (s$y_samples - y)
  return(list(
    id = i,
    y = y,
    y_bar = mean(s$y_samples),
    mse = mean(error_samples^2),
    mae = mean(abs(error_samples)),
    crps = crps(s$y_samples, y),
    t = get_time(fit)
  ))
}
