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
  arma::mat tmp;
  
  for (int r=0; r<s_many.n_rows; ++r) {
    tmp = s_few.row(r);
    s_many.row(r) = tmp * K;
  }
  
  /*
  arma::rowvec tmp;
  for (int i=0; i<s_many.n_rows; ++i) {
    for (int j=0; j<s_many.n_cols; ++j) {
      tmp = s_few.tube(i, j);
      s_many.tube(i, j) = tmp * K;     
    }
  }
  */
  
  return;
}