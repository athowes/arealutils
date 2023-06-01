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
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  
  nll -= dbinom_robust(y, m, beta_0, true).sum();
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
