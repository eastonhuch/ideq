#include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Samples from a multivariate normal distribution
//'
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

  // Calculate square root matrix of Sigma
  arma::mat Sigma_sqrt;
  try {
    Sigma_sqrt = arma::sqrtmat_sympd(Sigma);
  }
  catch (std::runtime_error e)
  {
    Rcout << "Failed to calculate sqrt(Sigma) using sqrtmat_sympd" << std::endl;
    Rcout << "Forcing symmetry using (Sigma + Sigma')/2 and trying again" << std::endl;
    arma::mat Sigma_sym = ( Sigma + Sigma.t() ) / 2;
    Sigma_sqrt = arma::sqrtmat_sympd(Sigma_sym);
  }

  arma::colvec x = mean + Sigma_sqrt * z;
  return x;
};

double rigamma(const double a, const double scl) {
  return (1 / R::rgamma(a, 1/scl));
}
