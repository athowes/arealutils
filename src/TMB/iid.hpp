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

  DATA_INTEGER(left_out); // Should data be left out
  DATA_IVECTOR(ii);       // Index of missing observations (zero-indexed)
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(u); // Spatial effects
  PARAMETER(log_sigma_u); // Log standard deviation of spatial effects
  
  // Transformed parameters block
  Type sigma_u(exp(log_sigma_u));
  vector<Type> eta(beta_0 + u);
  vector<Type> rho(invlogit(eta));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_u, Type(0), Type(2.5), true) + log_sigma_u; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  nll -= dnorm(u, Type(0), sigma_u, true).sum();
  
  vector<Type> log_lik(dbinom_robust(y, m, eta, true));
  
  // ADREPORT before zeroing some of the log_lik
  ADREPORT(log_lik);
  
  if(left_out > 0) {
    for (int j = 0; j < ii.size(); j++) {
      log_lik[ii[j]] = Type(0);
    }
  }
  
  nll -= log_lik.sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
