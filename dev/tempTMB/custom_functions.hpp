/// @file custom_functions.hpp

#ifndef custom_functions_hpp
#define custom_functions_hpp

template<class Type>
Type matrix_average(const matrix<Type> mat) {
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
    K(i, i) = matrix_average(block(kS, start_i, start_i, length_i, length_i));
    for (int j = i + 1; j < n; j++) {
      // Off-diagonal entries
      K(i, j) = matrix_average(block(kS, start_i, start_index(j), length_i, sample_lengths(j)));
      K(j, i) = K(i, j);
    }
  }
  K(n - 1, n - 1) = matrix_average(block(kS, start_index(n - 1), start_index(n - 1), sample_lengths(n - 1), sample_lengths(n - 1)));
  
  return K;
}

#endif
