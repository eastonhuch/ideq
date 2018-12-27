#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::colvec mvnorm(const arma::colvec & mean, const arma::mat & Sigma);

double ldmvnorm(const arma::mat & x, const arma::mat & mu, const arma::mat & Sigma);

double ldiwishart(const arma::cube & x, const double df,
                  const arma::cube & scale);

double rigamma(const double a, const double scl);

#endif
