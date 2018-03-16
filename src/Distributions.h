#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::colvec mvnorm(arma::colvec mean, arma::mat Sigma);

double rigamma(double a, double scl);

#endif
