% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmb.R
\name{ik_tmb}
\alias{ik_tmb}
\title{Fit Bayesian Integrated MVN Small Area Estimation model using \code{TMB}.}
\usage{
ik_tmb(sf, its = 1000, L = 10, type = "hexagonal", ii = NULL, ...)
}
\arguments{
\item{sf}{A simple features object with some geometry.}

\item{its}{Number of iterations in outer loop optimisation, passed to \code{nlminb}.}

\item{L}{The number of Monte Carlo samples to draw from each area.}

\item{type}{The \code{type} argument of \code{sf::st_sample}, defaults to \code{"hexagonal"}}

\item{ii}{The (zero-indexed) indices of the observations held-out.}

\item{...}{Additional arguments to \code{kernel}.}
}
\description{
Random effects have a multivariate Gaussian distribution with covariance
matrix calculated using \code{\link{integrated_covariance}}. Kernel hyperparameters
are given a prior and learnt.
}
\examples{
ik_tmb(mw, nsim_warm = 0, nsim_iter = 100, cores = 2)
}
