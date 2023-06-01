/// @file iid.hpp

#ifndef iid_hpp
#define iid_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type iid(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(phi); // Spatial effects
  PARAMETER(log_sigma_phi); // Log standard deviation of spatial effects
  
  // Transformed parameters block
  Type sigma_phi(exp(log_sigma_phi));
  vector<Type> eta(beta_0 + sigma_phi * phi);
  vector<Type> rho(invlogit(eta));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_phi, Type(0), Type(2.5), true) + log_sigma_phi; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  nll -= dnorm(phi, Type(0), sigma_phi, true).sum();
  
  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
