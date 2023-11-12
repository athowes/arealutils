// integrated.stan: Fully Bayesian integrated kernel plus CV

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
  
  // Inverse Gamma prior
  real<lower=0> a;
  real<lower=0> b;
  
  int sample_lengths[n]; // Number of Monte Carlo samples in each area
  int<lower=0> total_samples; // sum(sample_lengths)
  int start_index[n]; // Start indicies for each group of samples
  matrix[total_samples, total_samples] S; // Distances between all points (could be sparser!)
}

parameters {
  vector<lower=0>[n_mis] y_mis; // Vector of missing responses
  real beta_0; // Intercept
  vector[n] u; // Spatial effects
  real<lower=0> sigma_u; // Standard deviation of spatial effects
  real<lower=0> l; // Kernel lengthscale
}

transformed parameters {
  vector[n] eta = beta_0 + sigma_u * u;
  
  vector[n] y;
  y[ii_obs] = y_obs;
  y[ii_mis] = y_mis;
}

model {
  matrix[n, n] K = cov_sample_average(S, l, n, start_index, sample_lengths, total_samples);
  l ~ inv_gamma(a, b);
  sigma_u ~ normal(0, 2.5); // Weakly informative prior
  beta_0 ~ normal(-2, 1);
  u ~ multi_normal(mu, K);
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
