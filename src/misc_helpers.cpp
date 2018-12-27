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

