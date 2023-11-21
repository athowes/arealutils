/// @file besag.hpp

#ifndef besag_hpp
#define besag_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type besag(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR
  DATA_SCALAR(Qrank);
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(u); // Spatial effects
  PARAMETER(log_sigma_u); // Log standard deviation of spatial effects
  
  // Transformed parameters block
  Type sigma_u(exp(log_sigma_u));
  vector<Type> eta(beta_0 + u);
  vector<Type> rho(invlogit(eta));
  
  // Initialise negative log-likelihood
  Type nll;
  nll = Type(0.0);
  
  // Likelihood from priors
  nll -= dnorm(sigma_u, Type(0), Type(2.5), true) + log_sigma_u; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  
  // Besag
  nll -=  Qrank * 0.5 * log_sigma_u - 0.5 * sigma_u * (u * (Q * u)).sum();
  nll -= dnorm(u.sum(), Type(0.0), Type(0.001) * n, true); // Soft sum-to-zero constraint
  
  vector<Type> log_lik(dbinom_robust(y, m, eta, true));
  nll -= log_lik.sum();
  
  ADREPORT(log_lik);
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
