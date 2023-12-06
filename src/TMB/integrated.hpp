/// @file centroid.hpp

#ifndef integrated_hpp
#define integrated_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "../custom_functions.hpp"

template<class Type>
Type integrated(objective_function<Type>* obj) {
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  DATA_IVECTOR(ii_mis); // Index of missing observations (zero-indexed)
  
  // Inverse Gamma prior
  DATA_SCALAR(a);
  DATA_SCALAR(b);
  
  DATA_IVECTOR(sample_lengths); // Number of Monte Carlo samples in each area
  DATA_INTEGER(total_samples); // sum(sample_lengths)
  DATA_IVECTOR(start_index); // Start indices for each group of samples
  
  DATA_MATRIX(S); // Distances between all points (could be sparser!)
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(u); // Spatial effects
  PARAMETER(log_sigma_u); // Log standard deviation of spatial effects
  PARAMETER(log_l); // Log kernel lengthscale
  
  // Transformed parameters block
  Type sigma_u(exp(log_sigma_u));
  Type l(exp(log_l));
  vector<Type> eta(beta_0 + u);
  vector<Type> rho(invlogit(eta));
  matrix<Type> K(cov_sample_average(S, l, n, start_index, sample_lengths, total_samples));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_u, Type(0), Type(2.5), true) + log_sigma_u; // Change of variables
  nll -= dnorm(beta_0, Type(-2), Type(1), true); // NB: true puts the likelihood on the log-scale
  
  nll -= a * log(b) - lgamma(a) - (a + 1) * log(l) - b / l; // Inverse gamma
  nll -= log_l; // Change of variables
  
  using namespace density;
  K *= sigma_u;
  nll += MVNORM(K)(u); // On the negative log-scale already
  
  vector<Type> log_lik(dbinom_robust(y, m, eta, true));
  
  // ADREPORT before zeroing some of the log_lik
  ADREPORT(log_lik);
  
  if(ii_mis.size() > 0) {
    for (int i = 0; i < ii_mis.size(); i++) {
      log_lik[ii_mis[i]] = Type(0);
    }
  }
  
  nll -= log_lik.sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates
  
  return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
