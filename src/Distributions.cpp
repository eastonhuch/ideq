  // [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec mvnorm(arma::colvec mean, arma::mat Sigma) {
  int n = mean.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; i++) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }
  arma::colvec x = mean + arma::chol(Sigma).t() * z;
  return x;
};

double rigamma(double a, double scl) {
  return (1 / R::rgamma(a, 1 / scl));
}
