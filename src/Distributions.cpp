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

  // Calculate square root matrix of Sigma
  arma::mat Sigma_sqrt;
  try {
    Sigma_sqrt = arma::sqrtmat_sympd(Sigma);
  }
  catch (std::runtime_error e)
  {
    Rcout << "Failed to calculate sqrt(Sigma) using sqrtmat_sympd" << std::endl;
    Rcout << "Forcing symmetry using (Sigma + Sigma')/2 and trying again" << std::endl;
    Rcout << Sigma << std::endl;
    arma::mat Sigma_sym = ( Sigma + Sigma.t() ) / 2;
    Rcout << Sigma_sym << std::endl;
    arma::mat Sigma_sqrt = arma::sqrtmat_sympd(Sigma_sym);
    Rcout << Sigma_sqrt << std::endl;
  }

  arma::colvec x = mean + Sigma_sqrt * z;
  return x;
};

// Should this be the scale or rate parameterization?
double rigamma(const double a, const double scl) {
  return (1 / R::rgamma(a, 1/scl));
}
