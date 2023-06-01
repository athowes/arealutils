#include <TMB.hpp>

template<class Type>
Type mean(const matrix<Type> mat) {
  int num_elements = mat.rows() * mat.cols();
  Type sum = 0;
  
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      sum += mat(i, j);
    }
  }
  
  return sum / num_elements;
}

template<class Type>
matrix<Type> block(const matrix<Type> mat, int start_row, int start_col, int row_length, int col_length) {
  matrix<Type> submat(row_length, col_length);
  
  for (int i = 0; i < row_length; i++) {
    for (int j = 0; j < col_length; j++) {
      submat(i, j) = mat(start_row + i, start_col + j);
    }
  }
  
  return submat;
}

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

template<class Type>
matrix<Type> cov_sample_average(matrix<Type> S, Type l, int n, vector<int> start_index, vector<int> sample_lengths, int total_samples) {
  matrix<Type> K(n, n);
  matrix<Type> kS = cov_matern32(S, l);
  
  for (int i = 0; i < n - 1; i++) {
    // Diagonal entries (apart from the last)
    int start_i = start_index(i);
    int length_i = sample_lengths(i);
    K(i, i) = mean(block(kS, start_i, start_i, length_i, length_i));
    for (int j = i + 1; j < n; j++) {
      // Off-diagonal entries
      K(i, j) = mean(block(kS, start_i, start_index(j), length_i, sample_lengths(j)));
      K(j, i) = K(i, j);
    }
  }
  K(n - 1, n - 1) = mean(block(kS, start_index(n - 1), start_index(n - 1), sample_lengths(n - 1), sample_lengths(n - 1)));
  
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
  
  vector<int> sample_lengths(n); // Number of Monte Carlo samples in each area
  int total_samples; // sum(sample_lengths)
  vector<int> start_index(n); // Start indices for each group of samples
  matrix<Type> S(total_samples, total_samples); // Distances between all points (could be sparser!)
  
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
  matrix<Type> K(cov_sample_average(S, l, n, start_index, sample_lengths, total_samples));
  
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
