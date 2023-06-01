#include <TMB.hpp>

template<class Type>
matrix<Type> cov_matern32(matrix<Type> D, Type l) {
  int n = D.rows();
  matrix<Type> K(n, n);
  Type norm_K;
  Type sqrt3 = sqrt(3.0);
  
  for (int i = 0; i < n - 1; i++) {
    // Diagonal entries (apart from the last)
    K(i, i) = 1;
    for (int j = i + 1; j < n; j++) {
      // Off-diagonal entries
      norm_K = D(i, j) / l;
      K(i, j) = (1 + sqrt3 * norm_K) * exp(-sqrt3 * norm_K); // Fill lower triangle
      K(j, i) = K(i, j); // Fill upper triangle
    }
  }
  K(n - 1, n - 1) = 1;
  
  return K;
}

template <class Type>
Type objective_function<Type>::operator()()
{
  // Data block
  DATA_INTEGER(n); // Number of regions
  DATA_VECTOR(y); // Vector of responses
  DATA_VECTOR(m); // Vector of sample sizes
  DATA_MATRIX(D); // Distance between centroids
  
  // Inverse Gamma prior
  DATA_SCALAR(a);
  DATA_SCALAR(b);
  
  // Parameter block
  PARAMETER(beta_0); // Intercept
  PARAMETER_VECTOR(phi); // Spatial effects
  PARAMETER(log_sigma_phi); // Log standard deviation of spatial effects
  PARAMETER(log_l); // Log kernel lengthscale
  
  // Transformed parameters block
  Type sigma_phi(exp(log_sigma_phi));
  Type l(exp(log_l));
  vector<Type> eta(beta_0 + sigma_phi * phi);
  vector<Type> rho(invlogit(eta));
  matrix<Type> K(cov_matern32(D, l));
  
  // Model
  Type nll;
  nll = Type(0.0);
  
  nll -= dnorm(sigma_phi, Type(0), Type(100), true) + log_sigma_phi; // Approximating the uniform prior
  nll -= dnorm(beta_0, Type(-2), Type(5), true); // NB: true puts the likelihood on the log-scale
  
  nll -= a * log(b) - lgamma(a) - (a + 1) * log(l) - b / l; // Inverse gamma
  nll -= log_l; // Change of variables
  
  using namespace density;
  nll += MVNORM(K)(phi); // On the negative log-scale already

  nll -= dbinom_robust(y, m, eta, true).sum();
  
  ADREPORT(rho); // Would like to see posterior prevalence estimates

  return(nll);
}
