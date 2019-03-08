#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>

using namespace Rcpp;

void makeSymmetric(arma::mat & X) {
  X = (X + X.t())/2;
  return;
}

arma::mat forceInv(arma::mat X) {
  Rcout << "Forcing symmetry and trying again" << std::endl;
  makeSymmetric(X);
  return arma::inv_sympd(X);
}

arma::mat forceSqrtMat(arma::mat X) {
  Rcout << "Forcing symmetry and trying again" << std::endl;
  makeSymmetric(X);
  return arma::sqrtmat_sympd(X);
}