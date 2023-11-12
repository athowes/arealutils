/// @file bym2.hpp

#ifndef bym2_hpp
#define bym2_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type bym2(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  
  DATA_SPARSE_MATRIX(Q); // Structure matrix for ICAR
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(u); // Spatial effects
  PARAMETER_VECTOR(w); // Spatial component of spatial effects
  PARAMETER(logit_phi); // 
  PARAMETER(log_sigma_u); // Log standard deviation of spatial effects
  
  // Transformed parameters block
  Type sigma_u(exp(log_sigma_u));
  vector<Type> eta(beta_0 + sigma_u * u);
  vector<Type> rho(invlogit(eta));
  Type phi(invlogit(logit_phi));
  
  // Initialise negative log-likelihood
  Type nll;
  nll = Type(0.0);
  
  // Likelihood from priors
  nll -= dnorm(sigma_u, Type(0), Type(2.5), true) + log_sigma_u; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  
  nll -= log(phi) +  log(1 - phi);  // Change of variables: logit_phi -> phi
  nll -= dbeta(phi, Type(0.5), Type(0.5), true);
  
  // BYM2
  // Constant terms omitted: -0.5 * (n + rank(Q)) * log(2 * M_PI) + 0.5 * log|Q|
  nll -= -0.5 * n * (2 * log(sigma_u) + log(1 - phi));  // Normalising constant
  nll -= -0.5 / (sigma_u * sigma_u * (1 - phi)) * (u * u).sum();
  nll -= sqrt(phi) / (sigma_u * (1 - phi)) * (u * w).sum();
  nll -= -0.5 * (w * (Q * w)).sum();
  nll -= -0.5 * phi / (1 - phi) * (w * w).sum();
  
  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
