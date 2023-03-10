% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_marginal.R
\name{sample_marginal}
\alias{sample_marginal}
\alias{sample_marginal.inla}
\alias{sample_marginal.stanfit}
\title{Marginal samples.}
\usage{
sample_marginal(fit, ...)

\method{sample_marginal}{inla}(fit, i, n_obs, S = 1000, ...)

\method{sample_marginal}{stanfit}(fit, i, n_obs, S = 1000, ...)
}
\arguments{
\item{fit}{Fitted model.}

\item{...}{Additional arguments passed to \code{sample_marginal}.}

\item{i}{The index of the marginal to sample from.}

\item{n_obs}{The sample size.}

\item{S}{The number of draws.}
}
\description{
Marginal samples.
}
