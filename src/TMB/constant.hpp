/// @file constant.hpp

#ifndef constant_hpp
#define constant_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type constant(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes

  DATA_INTEGER(left_out); // Should data be left out
  DATA_IVECTOR(ii);       // Index of missing observations (zero-indexed)
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale

  vector<Type> log_lik(dbinom_robust(y, m, beta_0, true));
  
  // ADREPORT before zeroing some of the log_lik
  ADREPORT(log_lik);
  
  if(left_out > 0) {
    for (int j = 0; j < ii.size(); j++) {
      log_lik[ii[j]] = Type(0);
    }
  }
  
  nll -= log_lik.sum();

  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
