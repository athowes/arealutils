#include <TMB.hpp>
#include "custom_functions.hpp"

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  // Inverse Gamma prior
  DATA_SCALAR(a);
  DATA_SCALAR(b);
  
  DATA_IVECTOR(sample_lengths); // Number of Monte Carlo samples in each area
  DATA_INTEGER(total_samples); // sum(sample_lengths)
  DATA_IVECTOR(start_index); // Start indices for each group of samples
  
  DATA_MATRIX(S); // Distances between all points (could be sparser!)
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(phi); // Spatial effects
  PARAMETER(log_sigma_phi); // Log standard deviation of spatial effects
  PARAMETER(log_l); // Log kernel lengthscale
  
  // Transformed parameters block
  Type sigma_phi(exp(log_sigma_phi));
  Type l(exp(log_l));
  vector<Type> eta(beta_0 + sigma_phi * phi);
  vector<Type> rho(invlogit(eta));
  matrix<Type> K(cov_sample_average(S, l, n, start_index, sample_lengths, total_samples));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_phi, Type(0), Type(2.5), true) + log_sigma_phi; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  
  nll -= a * log(b) - lgamma(a) - (a + 1) * log(l) - b / l; // Inverse gamma
  nll -= log_l; // Change of variables

  nll += density::MVNORM(K)(phi); // On the negative log-scale already
  
  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}
