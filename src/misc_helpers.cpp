#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void makeSymmetric(arma::mat & X) {
  X = (X + X.t())/2;
  return;
}

arma::mat forceInv(arma::mat X) {
  arma::mat X_inv;
  
  Rcout << "Forcing symmetry and trying again" << std::endl;
  makeSymmetric(X);
  try {
    X_inv = arma::inv_sympd(X);
  } 
  catch (std::runtime_error e) {
    Rcout << "Inversion unsuccessful" << std::endl; 
    Rcout << "Adding 1 to diagonal and trying again" << std::endl;
    X.diag() += 1;
    X_inv = arma::inv_sympd(X);
    Rcout << "Inversion successful but sampling has not converged" << std::endl;
  }
  
  return X_inv;
}

arma::mat forceSqrtMat(arma::mat X) {
  arma::mat X_sqrt;

  Rcout << "Forcing symmetry and trying again" << std::endl;
  makeSymmetric(X);
  try {
    X_sqrt = arma::sqrtmat_sympd(X);
    Rcout << "sqrt(X) successful" << std::endl;
  }
  catch (std::runtime_error e)
  {
    Rcout << "sqrt(X) unsuccessful" << std::endl;
    Rcout << "Adding 1 to diagonal and trying again" << std::endl;
    X.diag() += 1;
    X_sqrt = arma::sqrtmat_sympd(X);
    Rcout << "sqrt(X) successful but sampling has not converged" << std::endl;
  }
  
  return X_sqrt;
}