  // [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

//' Samples from a multivariate normal distribution
//'
//' @keyword Normal, Multivariate, Sample, Random
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
arma::colvec mvnorm(const arma::colvec & mean, const arma::mat & Sigma) {
  int n = mean.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; ++i) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }
  arma::colvec x = mean + arma::chol(Sigma).t() * z;
  return x;
};

double rigamma(const double a, const double scl) {
  return (1 / R::rgamma(a, 1 / scl));
}
