/// @file mvn_covariance.hpp

#ifndef mvn_covariance_hpp
#define mvn_covariance_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type mvn_covariance(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  DATA_MATRIX(Sigma); // Dense covariance matrix
  
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
  
  using namespace density;
  Sigma *= sigma_u;
  nll += MVNORM(Sigma)(u); // On the negative log-scale already

  vector<Type> log_lik(dbinom_robust(y, m, eta, true));
  nll -= log_lik.sum();
  
  ADREPORT(log_lik);
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
