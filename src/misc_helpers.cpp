#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void make_symmetric(arma::mat & X) {
  const int n = X.n_rows;
  double tmp = 0;
  for (int i=1; i<n; ++i) {
    for (int j=0; j<i; ++j) {
      tmp = (X.at(i, j) + X.at(j, i)) / 2;
      X.at(i, j) = tmp;
      X.at(j, i) = tmp;
    }
  }
  return;
};

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K) {
  if (s_many.n_rows != s_few.n_rows || s_many.n_cols != s_few.n_cols) {
    throw std::invalid_argument("s_many and s_few must have same number of rows and columns");
  }
  arma::colvec tmp;
  
  for (int r=0; r<s_few.n_rows; ++r) {
    for (int c=0; c<s_few.n_cols; ++c) {
      tmp = s_few.tube(r, c);
      s_many.tube(r, c) = K * tmp;
    }
  }
  
  return;
}