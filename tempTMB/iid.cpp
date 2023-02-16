#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(phi); // Spatial effects
  PARAMETER(sigma_phi); // Standard deviation of spatial effects
  
  // Transformed parameters block
  vector<Type> eta(beta_0 + sigma_phi * phi);
  vector<Type> rho(invlogit(eta));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_phi, Type(0), Type(100), true); // Approximating the uniform prior
  nll -= dnorm(beta_0, Type(-2), Type(5), true); // NB: true puts the likelihood on the log-scale
  nll -= dnorm(phi, Type(0), sigma_phi, true).sum();
  
  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}
