#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR

  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(phi); // Spatial effects
  PARAMETER_VECTOR(phi_u); // Spatial component of spatial effects
  PARAMETER(logit_pi); // 
  PARAMETER(sigma_phi); // Standard deviation of spatial effects
  
  // Transformed parameters block
  vector<Type> eta(beta_0 + sigma_phi * phi);
  vector<Type> rho(invlogit(eta));
  Type pi(invlogit(logit_pi));
  
  // Initialise negative log-likelihood
  Type nll;
  nll = Type(0.0);
  
  // Likelihood from priors
  nll -= dnorm(sigma_phi, Type(0), Type(100), true); // Approximating the uniform prior
  nll -= dnorm(beta_0, Type(-2), Type(5), true); // NB: true puts the likelihood on the log-scale
  
  nll -= log(pi) +  log(1 - pi);  // Change of variables: logit_pi -> pi
  nll -= dbeta(pi, Type(0.5), Type(0.5), true);
  
  // BYM2
  // Constant terms omitted: -0.5 * (n + rank(Q)) * log(2 * M_PI) + 0.5 * log|Q|
  nll -= -0.5 * n * (2 * log(sigma_phi) + log(1 - pi));  // Normalising constant
  nll -= -0.5 / (sigma_phi * sigma_phi * (1 - pi)) * (phi * phi).sum();
  nll -= sqrt(pi) / (sigma_phi * (1 - pi)) * (phi * phi_u).sum();
  nll -= -0.5 * (phi_u * (Q * phi_u)).sum();
  nll -= -0.5 * pi / (1 - pi) * (phi_u * phi_u).sum();
  
  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}
