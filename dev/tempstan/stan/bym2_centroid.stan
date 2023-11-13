// bym2_centroid.stan: Fully Bayesian centroid kernel with IID noise plus CV

functions{
#include /include/custom_functions.stan
}

data {
  int<lower=0> n_obs; // Number of observed regions
  int<lower=0> n_mis; // Number of missing regions
  
  int<lower = 1, upper = n_obs + n_mis> ii_obs[n_obs];
  int<lower = 1, upper = n_obs + n_mis> ii_mis[n_mis];
  
  int<lower=0> n; // Number of regions n_obs + n_mis
  vector[n_obs] y_obs; // Vector of observed responses
  vector[n] m; // Vector of sample sizes
  vector[n] mu; // Prior mean vector
  matrix[n, n] D; // Distances between centroids
}

parameters {
  vector<lower=0>[n_mis] y_mis; // Vector of missing responses
  real beta_0; // Intercept
  vector[n] w; // Structured spatial effects
  vector[n] v; // Unstructured spatial effects
  real<lower=0, upper=1> phi; // Proportion unstructured vs. structured variance
  real<lower=0> sigma_u; // Standard deviation of spatial effects
  real<lower=0> l; // Kernel lengthscale
}

transformed parameters {
  vector[n] u = sqrt(1 - phi) * v + sqrt(phi) * w; // Spatial effects
  vector[n] eta = beta_0 + sigma_u * u;
  
  vector[n] y;
  y[ii_obs] = y_obs;
  y[ii_mis] = y_mis;
}

model {
  matrix[n, n] K = cov_matern32(D, l);
  // I could do this?
  // matrix[n, n] L = cholesky_decompose(K);
  // y ~ multi_normal_cholesky(mu, L);
  l ~ gamma(1, 1);
  sigma_u ~ normal(0, 2.5); // Weakly informative prior
  beta_0 ~ normal(-2, 1);
  phi ~ beta(0.5, 0.5);
  v ~ normal(0, 1);
  w ~ multi_normal(mu, K);
  for(i in 1:n) {
    y[i] ~ xbinomial_logit(m[i], eta[i]); 
  }
}

generated quantities {
  real tau_u = 1 / sigma_u^2; // Precision of spatial effects
  vector[n] rho = inv_logit(beta_0 + sigma_u * u);
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = xbinomial_logit_lpdf(y[i] | m[i], eta[i]);
  }
}