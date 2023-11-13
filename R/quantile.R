#' Compute an the quantile of an observation within a vector of samples.
#'
#' @param samples A collection of draws from a probability distribution.
#' @param true_value The underlying true value that is being compared to.
#' @return The quantile.
#' @export
quantile <- function(samples, true_value){
  min(which(sort(samples) >= true_value)) / length(samples)
}
