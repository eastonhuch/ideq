#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::colvec mvnorm(const arma::colvec & mean, const arma::mat & Sigma);

double rigamma(const double a, const double scl);

#endif
